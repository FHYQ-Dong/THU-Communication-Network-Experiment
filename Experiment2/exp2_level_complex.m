% 通信与网络 实验2 基带传输
% 电平传输 复电平信道
clearvars;
close all;
clc;
rng(2025);  %固定种子

% 参数
% 信噪比范围 E_s/\sigma^2
SNR_dB = 0:1:25;
SNR = 10.^(SNR_dB/10);

% 符号数
M = [4;16];

% 重复实验次数，可修改
N = 5e3;
if mod(N,4)~=0
    error('Invalid bit sequence length');
end

% 记录BER和SER
BER = zeros(length(M),length(SNR));
SER = zeros(length(M),length(SNR));

% 记录两个信噪比下的接收信号
SNR_dB_rec = [15,20];

for cnt1 = 1:length(M)
    M_cur = M(cnt1);
    % 随机生成比特序列（data_bit）
    data_bit = randi([0 1], 1, N);

    % 符号映射，转化为调制后的电平序列（data_mod）
    if M_cur == 4
        [data_mod, mod_set] = Gray_mapping_PSK_complex(data_bit, M_cur);
    elseif M_cur == 16
        [data_mod, mod_set] = Gray_mapping_QAM_complex(data_bit, M_cur);
    end

    % 发送端平均电平能量
    Es = mean(abs(data_mod).^2);
    % 理论平均电平能量
    if M_cur == 4
        Es_theory = 1;
    elseif M_cur == 16
        Es_theory = 2*(M_cur-1)/3;
    end

    for cnt2 = 1:length(SNR)
        % 加性噪声
        noise_power = Es_theory ./ SNR(cnt2);
        
        noise = sqrt(noise_power./2) .* (randn(1, N/log2(M_cur)) + 1j .* randn(1, N/log2(M_cur))); 

        % 加性高斯噪声信道，获得接收信号（data_rx）
        data_rx = data_mod + noise;
        if M_cur==16 && ismember(SNR_dB(cnt2),SNR_dB_rec)
            figure;
            scatter(real(data_rx), imag(data_rx), 'b.');
            hold on;
            [mod_set_r, mod_set_i] = meshgrid(real(mod_set), imag(mod_set));
            scatter(mod_set_r, mod_set_i, 'ro');
            box on; grid on;
            title(['SNR=',num2str(SNR_dB(cnt2)),'dB']);
            hold off;
        end

        % 解调，根据电平序列获得比特序列（data_bit_demod）
        if M_cur == 4
            data_bit_demod = Inverse_Gray_mapping_PSK_complex(data_rx, M_cur, mod_set);
        elseif M_cur == 16
            data_bit_demod = Inverse_Gray_mapping_QAM_complex(data_rx, M_cur, mod_set);
        end

        % 差错比特
        bit_err = xor(data_bit, data_bit_demod);
        % 计算差错比特数
        BER(cnt1,cnt2) = sum(bit_err);
        % 计算差错符号数
        bit_err = reshape(bit_err, log2(M_cur), []);
        SER(cnt1,cnt2) = sum(any(bit_err,1));
    end
end
% 误码率和误比特率归一化
BER = BER./N;
SER = SER ./ (N./log2(M)); 

% 理论值
ser_theory = [2.*myqfunc(sin(pi/M(1)).*sqrt(SNR.*2)); 4.*(1-1/sqrt(M(2))).*myqfunc(sqrt(3/2/(M(2)-1).*SNR))];
ber_theory = ser_theory./log2(M);

% 绘图
figure;
disp("size of SNR_dB"); disp(size(SNR_dB));
disp("size of BER"); disp(size(BER));
disp("size of ber_theory"); disp(size(ber_theory));
disp(BER(2,:));
semilogy(SNR_dB, BER(1,:),'bo', ...
    SNR_dB, ber_theory(1,:),'b', ...
    SNR_dB, BER(2,:),'rx', ...
    SNR_dB, ber_theory(2,:),'r');
xlabel('信噪比E_s/\sigma^2 (dB)');
ylabel('误比特率');
set(gca, 'ylim', [1e-4,1e0]);
legend('M=4, simulation', 'M=4, theory', ...
    'M=16, simulation', 'M=16, theory', ...
    'Location', 'southwest');
title('复电平信道 BER-SNR');
box on;
grid on;

figure;
semilogy(SNR_dB, SER(1,:),'bo', ...
    SNR_dB, ser_theory(1,:),'b', ...
    SNR_dB, SER(2,:),'rx', ...
    SNR_dB, ser_theory(2,:),'r');
xlabel('信噪比E_s/\sigma^2 (dB)');
ylabel('误符号率');
set(gca, 'ylim', [1e-4,1e0]);
legend('M=4, simulation', 'M=4, theory', ...
    'M=16, simulation', 'M=16, theory', ...
    'Location', 'southwest');
title('复电平信道 SER-SNR');
box on;
grid on;


function [data_mod, mod_set] = Gray_mapping_QAM_complex(data_bit, M)
    % QAM Gray 编码为电平序列（A = 1）
    %
    % 输入
    %   data_bit    待编码 bit 串（向量）
    %   M           星座图电平数（标量）
    %
    % 输出
    %   data_mod    复电平序列（行向量）
    %   mod_set     电平符号集合（复行向量）

    bit_per_sym = log2(M);
    data_bit = data_bit(:)';  % 转为行向量
    bit_mat = reshape(data_bit, bit_per_sym, []);  % 每列一个 Gray 码
    bits_I = bit_mat(1:bit_per_sym/2, :);  % 前 bit_per_sym 个 bit 为 I 路
    bits_Q = bit_mat(bit_per_sym/2+1:end, :);  % 后 bit_per_sym 个 bit 为 Q 路
    sym_I = 2 * (2.^((bit_per_sym/2-1):-1:0) * double(mod(cumsum(bits_I, 1), 2))) - sqrt(M) + 1;  % PAM 电平映射
    sym_Q = 2 * (2.^((bit_per_sym/2-1):-1:0) * double(mod(cumsum(bits_Q, 1), 2))) - sqrt(M) + 1;  % PAM 电平映射
    data_mod = complex(sym_I, sym_Q);  % 合并 I,Q 路
    mod_set = complex(-(sqrt(M) - 1):2:(sqrt(M) - 1), -(sqrt(M) - 1):2:(sqrt(M) - 1));  % QAM 电平集合
end


function [data_mod, mod_set] = Gray_mapping_PSK_complex(data_bit, M)
    % PSK Gray 编码为电平序列（A = 1）
    %
    % 输入
    %   data_bit    待编码 bit 串（向量）
    %   M           星座图电平数（标量）
    %
    % 输出
    %   data_mod    复电平序列（行向量）
    %   mod_set     电平符号集合（复行向量）

    bit_per_sym = log2(M);
    data_bit = data_bit(:).';  % 转为行向量
    data_mat = reshape(data_bit, bit_per_sym, []);  % 每列是一个 gray 码
    data_mod = exp(2j*pi/M .* (2.^((bit_per_sym-1):-1:0)) * double(mod(cumsum(data_mat, 1), 2)));  % PSK 电平映射
    mod_set = exp(2j*pi/M .* (0:M-1));  % PSK 电平集合
end


function data_bit_demod = Inverse_Gray_mapping_QAM_complex(data_rx, M, mod_set)
    % 硬判决并 QAM Gray 解码
    %
    % 输入
    %   data_rx           接收的电平序列（向量）
    %   M                 星座图电平数（标量）
    %   mod_set           电平符号集合（复向量）
    %
    % 输出
    %   data_bit_demod    解码的 bit 串（行向量）
    
    bit_per_sym = log2(M);
    I_rx = real(data_rx);  % 分开 I,Q 路
    Q_rx = imag(data_rx);  % 分开 I,Q 路
    I_mod_idx = dsearchn(real(mod_set(:)), I_rx(:));  % ML 准则找电平的索引，为一个列向量
    Q_mod_idx = dsearchn(imag(mod_set(:)), Q_rx(:));  % ML 准则找电平的索引，为一个列向量
    I_idx_bit_mat = zeros(bit_per_sym/2, length(I_rx));  % 索引的二进制表示，每列为一个索引
    Q_idx_bit_mat = zeros(bit_per_sym/2, length(Q_rx));  % 索引的二进制表示，每列为一个索引
    for k = 1:bit_per_sym/2 
        I_idx_bit_mat(k, :) = bitget(I_mod_idx, bit_per_sym/2-k+1).';
        Q_idx_bit_mat(k, :) = bitget(Q_mod_idx, bit_per_sym/2-k+1).';
    end
    I_gray_mat = I_idx_bit_mat;
    Q_gray_mat = Q_idx_bit_mat;
    I_gray_mat(2:end, :) = xor(I_idx_bit_mat(2:end, :), I_idx_bit_mat(1:end-1, :));  % 差分计算 Gray 码
    Q_gray_mat(2:end, :) = xor(Q_idx_bit_mat(2:end, :), Q_idx_bit_mat(1:end-1, :));  % 差分计算 Gray 码
    gray_mat = [I_gray_mat; Q_gray_mat];  % 合并 I,Q 路
    data_bit_demod = reshape(gray_mat, 1, []);  % 拉伸回 1 维
end


function data_bit_demod = Inverse_Gray_mapping_PSK_complex(data_rx, M, mod_set)
    % 硬判决并 PSK Gray 解码
    %
    % 输入
    %   data_rx           接收的电平序列（向量）
    %   M                 星座图电平数（标量）
    %   mod_set           电平符号集合（复向量）
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
