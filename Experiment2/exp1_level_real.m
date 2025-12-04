% 通信与网络 实验2 基带传输
% 电平传输 实电平信道
clearvars;
close all;
clc;
rng(2025);  %固定种子

% 参数
% 信噪比范围 E_s/\sigma^2
SNR_dB = -5:1:20;
SNR = 10.^(SNR_dB/10);
SNR_rec_id = [9, 13, 16]; % SNR = [2, 5, 10]

% 符号数
M = [2;4];

% 比特序列长度，可修改
N = 1e6;
if mod(N,2)~=0
    error('Invalid bit sequence length');
end

% 初始化向量，存储误符号率和误比特率
BER = zeros(length(M), length(SNR));
SER = zeros(length(M), length(SNR));

for cnt1 = 1:length(M)
    M_cur = M(cnt1);
    % 随机生成比特序列（data_bit）
    data_bit = randi([0 1], 1, N);

    % 符号映射，转化为调制后的电平序列（data_mod）
    [data_mod, mod_set] = Gray_mapping_real(data_bit, M_cur);

    % 发送端平均电平能量
    Es = mean(data_mod.^2, 'all');

    % 理论平均电平能量
    Es_theory = (M_cur^2-1)/3;

    for cnt2 = 1:length(SNR)
        % 加性噪声
        noise_power = Es_theory ./ SNR(cnt2);
        noise = sqrt(noise_power) .* randn(1, N/log2(M_cur));  % Ensure noise is a column vector

        % 加性高斯噪声信道，获得接收信号（data_rx）
        data_rx = data_mod + noise;

        % 解调，根据电平序列获得比特序列（data_bit_demod）
        data_bit_demod = Inverse_Gray_mapping(data_rx, M_cur, mod_set);

        % 差错比特
        bit_err = xor(data_bit, data_bit_demod);
        % 计算差错比特数
        BER(cnt1,cnt2) = sum(bit_err);
        % 计算差错符号数
        if M_cur == 2
            SER(cnt1,cnt2) = BER(cnt1,cnt2);
        elseif M_cur == 4
            bit_err = reshape(bit_err,2,[]);
            SER(cnt1,cnt2) = sum(bit_err(1,:) | bit_err(2,:));
        end
    end
end
% 误码率和误比特率归一化
BER = BER ./ N;
SER = SER ./ (N./log2(M)); 

% 理论值
SER_theory = (2.*M-2)./M .* myqfunc(sqrt(3.*SNR./(M.^2-1)));
BER_theory = SER_theory./log2(M);

% 绘图
figure;
semilogy(SNR_dB, BER(1,:),'bo', ...
    SNR_dB, BER_theory(1,:),'b', ...
    SNR_dB, BER(2,:),'rx', ...
    SNR_dB, BER_theory(2,:),'r');
xlabel('信噪比E_s/\sigma^2 (dB)');
ylabel('误比特率');
set(gca, 'ylim', [1e-4,1e0]);
legend('M=2, simulation', 'M=2, theory', ...
    'M=4, simulation', 'M=4, theory', ...
    'Location', 'southwest');
title('实电平信道 BER-SNR');
box on;
grid on;

figure;
semilogy(SNR_dB, SER(1,:),'bo', ...
    SNR_dB, SER_theory(1,:),'b', ...
    SNR_dB, SER(2,:),'rx', ...
    SNR_dB, SER_theory(2,:),'r');
xlabel('信噪比E_s/\sigma^2 (dB)');
ylabel('误符号率');
set(gca, 'ylim', [1e-4,1e0]);
legend('M=2, simulation', 'M=2, theory', ...
    'M=4, simulation', 'M=4, theory', ...
    'Location', 'southwest');
title('实电平信道 SER-SNR');
box on;
grid on;

% 打印信息
SNR_rec = SNR(SNR_rec_id);
SER_rec_M2 = SER(1, SNR_rec_id);
SER_rec_M4 = SER(2, SNR_rec_id);
BER_rec_M2 = BER(1, SNR_rec_id);
BER_rec_M4 = BER(2, SNR_rec_id);
fprintf('=====  仿真长度 N = %d  =====\n', N);
disp('-----  差错统计（M=2）  -----');
for k = 1:length(SNR_rec)
    fprintf('SNR = %.2f, 误符号数 = %d, 误bit数 = %d\n', SNR_rec(k), SER_rec_M2(k) .* (N./log2(M(1))), BER_rec_M2(k) .* N);
end
disp('-----  差错统计（M=4）  -----');
for k = 1:length(SNR_rec)
    fprintf('SNR = %.2f, 误符号数 = %d, 误bit数 = %d\n', SNR_rec(k), SER_rec_M4(k) .* (N./log2(M(2))), BER_rec_M4(k) .* N);
end
fprintf('\n');

function [data_mod, mod_set] = Gray_mapping_real(data_bit, M)
    % Gray 编码为电平序列（A = 1)
    %
    % 输入
    %   data_bit    待编码 bit 串（向量）
    %   M           星座图电平数（标量）
    %
    % 输出
    %   data_mod    电平序列（行向量）
    %   mod_set     电平符号集合（行向量）
    
    bit_per_sym = log2(M);
    data_bit = data_bit(:).';  % 转为行向量
    data_mat = reshape(data_bit, bit_per_sym, []);  % 每列是一个 gray 码
    data_mod = 2 * (2.^((bit_per_sym-1):-1:0)) * double(mod(cumsum(data_mat, 1), 2)) - M + 1;  % PAM 电平映射
    mod_set = 2*(0:(M-1)) - M + 1;  % PAM 电平集合
end


function data_bit_demod = Inverse_Gray_mapping(data_rx, M, mod_set)
    % 硬判决并 Gray 解码
    %
    % 输入
    %   data_rx           接收的电平序列（向量）
    %   M                 星座图电平数（标量）
    %   mod_set           电平符号集合（向量）
    %
    % 输出
    %   data_bit_demod    解码的 bit 串（行向量）

    bit_per_sym = log2(M);
    mod_idx = dsearchn(mod_set(:), data_rx(:)) - 1;  % ML 准则找电平的索引，为一个列向量
    idx_bit_mat = zeros(bit_per_sym, length(data_rx));  % 索引的二进制表示，每列为一个索引
    for k = 1:bit_per_sym
        idx_bit_mat(k, :) = bitget(mod_idx, bit_per_sym-k+1).';
    end
    gray_mat = idx_bit_mat;
    gray_mat(2:end, :) = xor(idx_bit_mat(2:end, :), idx_bit_mat(1:end-1, :));  % 差分计算 Gray 码
    data_bit_demod = reshape(gray_mat, 1, []);  % 拉伸回 1 维
end
