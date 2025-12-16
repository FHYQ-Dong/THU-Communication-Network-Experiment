% 通信与网络 实验3 载波传输
%4PAM
clearvars;
close all;
clc;


M = 4;
N = 1e5;        % 符号数
Es = 1;         % 符号能量
T = 0.01;       % 符号周期
A = sqrt(3*Es/(T*(M^2-1))); % 幅度
f_c = 500;      % 载波频率
delta_t = 1e-4; % 采样间隔
num = round(T/delta_t);
t = delta_t:delta_t:N*T;

% 成形脉冲
p_t = A*ones(1,num);

% 生成随机比特并映射到4PAM符号
bit_data = randi([0 1],1,2*N);
mod_data = my_gray_map_real_M4(bit_data);

% 基带成形信号
x0_t = zeros(1,N*num);
for k = 1:N
    x0_t((k-1)*num+1:k*num) = mod_data(k)*p_t;
end

% 调制为载波信号
x_t = x0_t .* (sqrt(2)*cos(2*pi*f_c*t));



%% ====================== 匹配滤波器 ======================
h1_t = sqrt(2/T) * cos(2*pi*f_c*(delta_t:delta_t:T)); % 方法1
h2_t = sqrt(2/T) * ones(1,num);                       % 方法2

% 最佳采样时刻 t=kT
t_sample = T:T:N*T;
idx_sample = round(t_sample/delta_t);

%% ====================== 仿真不同信噪比 ======================
EbN0_dB = 2:2:8;       % 信噪比范围 dB
EbN0 = 10.^(EbN0_dB/10);
Eb = Es / log2(M);   
N0_list = Eb ./ EbN0;
EsN0 = EbN0 * log2(M); 
EsN0_dB = 10.*log10(EsN0);

ser_theory = zeros(size(EbN0));
ber_theory = zeros(size(EbN0));
%% 不加噪声观察匹配滤波
y_t = x_t;

% % 匹配滤波输出
mf_out1 = conv_no_delay(y_t, fliplr(h1_t), delta_t);
% 方法2先下变频
y_base = y_t .* cos(2*pi*f_c*t);
mf_out2 = conv_no_delay(y_base, fliplr(h2_t), delta_t);     

% 前5个符号输出波形
figure;
subplot(2,1,1)
plot(t(1:round(5*T/delta_t)), mf_out1(1:round(5*T/delta_t)), 'LineWidth', 1);
hold on;
stem(t(idx_sample(1:5)), mf_out1(idx_sample(1:5)), 'r', 'filled');
xlabel('时间 (s)');
ylabel('幅度');
title('方法1：y(t)*余弦载波匹配滤波输出');
legend('匹配滤波输出（不加噪声）','最佳采样时刻','Location','best');
box on; grid on;

subplot(2,1,2)
plot(t(1:round(5*T/delta_t)), mf_out2(1:round(5*T/delta_t)), 'LineWidth', 1);
hold on;
stem(t(idx_sample(1:5)), mf_out2(idx_sample(1:5)), 'r', 'filled');
xlabel('时间 (s)');
ylabel('幅度');
title('方法2：y(t)先乘余弦后与方波匹配滤波输出');
legend('匹配滤波输出（不加噪声）','最佳采样时刻','Location','best');
box on; grid on;

%% 采样点偏移观察

offset_list = [0 1 2]; % 采样点偏离+-2，+-1，0

BER_sim1 = zeros(length(EbN0_dB),length(offset_list));   % 方法1仿真BER
BER_sim2 = zeros(length(EbN0_dB),length(offset_list));   % 方法2仿真BER
SER_sim1 = zeros(length(EbN0_dB),length(offset_list));   % 方法1仿真SER
SER_sim2 = zeros(length(EbN0_dB),length(offset_list));   % 方法2仿真SER

for idx = 1:length(N0_list)
    n0 = N0_list(idx);

    % AWGN噪声
    noise_power = n0/2/delta_t;
    noise = sqrt(noise_power) * randn(size(t));

    % 接收信号
    y_t = x_t + noise;

    % 匹配滤波
    mf1_out = conv_no_delay(y_t, fliplr(h1_t), delta_t);
    y_base = y_t .* cos(2*pi*f_c*t);
    mf2_out = conv_no_delay(y_base, fliplr(h2_t), delta_t);

    % 遍历采样偏差
    for j = 1:length(offset_list)
        offset_time = offset_list(j);

        idx_s = idx_sample + offset_time;

        new_ind = idx_s(1:length(idx_s)-1);

        y1 = mf1_out( new_ind);
        y2 = mf2_out( new_ind);

        x_hat1 = zeros(1,2*N-2);
        x_hat2 = zeros(1,2*N-2);
        for k = 1:N-1
            % 判决时归一化，除以 A*sqrt(T) 得到 {-3, -1, 1, 3}
            x_hat1((k-1)*2+1:k*2) = my_inverse_gray_map_real_M4(y1(k)/(A*sqrt(T)));
            x_hat2((k-1)*2+1:k*2) = my_inverse_gray_map_real_M4(y2(k)/(A*sqrt(T)));
        end

        % BER/SER 统计
        bit_err1 = xor(bit_data(1:2*N-2), x_hat1);
        bit_err2 = xor(bit_data(1:2*N-2), x_hat2);
        BER_sim1(idx,j) = sum(bit_err1)/length(bit_data);
        BER_sim2(idx,j) = sum(bit_err2)/length(bit_data);

        bit_err1_reshape = reshape(bit_err1,2,N-1);
        bit_err2_reshape = reshape(bit_err2,2,N-1);
        SER_sim1(idx,j) = sum(any(bit_err1_reshape, 1)) / (N-1);
        SER_sim2(idx,j) = sum(any(bit_err2_reshape, 1)) / (N-1);
    end
    ser_theory(idx)= (2*(M-1)/M) * qfunc(sqrt( (3/(M^2-1)) * 2 * EsN0(idx) ));
    ber_theory(idx)= ser_theory(idx) / log2(M);
end


%% ====================== 绘图 ======================
colors = lines(length(offset_list)); % 生成直观颜色
markers = {'o','s','d'};     % 不同偏差标记
figure;
hold on;
set(gca, 'YScale', 'log');
for j = 1:length(offset_list)
    semilogy(EbN0_dB, BER_sim1(:,j), '-', 'Color', colors(j,:), 'Marker', markers{j}, 'LineWidth',1.5);
end

for j = 1:length(offset_list)
    semilogy(EbN0_dB, BER_sim2(:,j), '--', 'Color', colors(j,:), 'Marker', markers{j}, 'LineWidth',1.5);
end

legend_entries = {};

for j = 1:length(offset_list)
    legend_entries{end+1} = sprintf('方法1 offset %d 个采样点', offset_list(j));
end
for j = 1:length(offset_list)
    legend_entries{end+1} = sprintf('方法2 offset %d 个采样点', offset_list(j));
end
semilogy(EbN0_dB, ber_theory, '-', 'Color', 'r',  'LineWidth',1.5);
legend_entries{end+1} = '理论值';

legend(legend_entries,'Location','best');
xlabel('E_b/N_0 (dB)');
ylabel('误比特率 (BER)');
title('4PAM误比特率性能曲线');
grid on; box on;

figure;hold on;
set(gca, 'YScale', 'log');
for j = 1:length(offset_list)
    semilogy(EsN0_dB, SER_sim1(:,j), '-', 'Color', colors(j,:), 'Marker', markers{j}, 'LineWidth',1.5);
end

for j = 1:length(offset_list)
    semilogy(EsN0_dB, SER_sim2(:,j), '--', 'Color', colors(j,:), 'Marker', markers{j}, 'LineWidth',1.5);
end


legend_entries = {};
for j = 1:length(offset_list)
    legend_entries{end+1} = sprintf('方法1 offset %d 个采样点', offset_list(j));
end
for j = 1:length(offset_list)
    legend_entries{end+1} = sprintf('方法2 offset %d 个采样点', offset_list(j));
end

semilogy(EsN0_dB, ser_theory, '-', 'Color', 'r',  'LineWidth',1.5);
legend_entries{end+1} = '理论值';

legend(legend_entries,'Location','best');
xlabel('E_s/N_0 (dB)');
ylabel('误码率 (BER)');
title('4PAM误码率性能曲线');

grid on; box on;

%% 发送信号功率谱
% 设置滑动窗口参数
window_length = 1024; % 窗口长度
overlap = 512; % 重叠长度
nfft = 1024; % FFT长度

% 初始化PSD估计
psdEstimate = zeros(nfft, 1);  % Include both positive and negative frequencies

% 滑动窗口处理
for start = 1:(window_length - overlap):(length(x_t) - window_length)
    % 提取窗口
    windowedSegment = x_t(start:start+window_length-1);
    
    % 窗函数处理（这里使用汉明窗）
    windowedSegment = windowedSegment .* hamming(window_length)';
    
    % 计算傅里叶变换
    Y = fft(windowedSegment, nfft);
    
    % 计算功率谱密度
    P2 = abs(Y/nfft).^2 / 4;
    
    % 将功率谱加入估计值（包括负频率）
    psdEstimate = psdEstimate + P2';
end

% 完成平均
psdEstimate = psdEstimate / ((length(x_t) - window_length) / (window_length - overlap));

% 计算频率轴
f = (-nfft/2:nfft/2-1) / (nfft * delta_t);  

%功率谱理论值
S_X = (A^2*T*(M^2-1)/6) * ((sinc((f-f_c)*T)).^2 + (sinc((f+f_c)*T)).^2);

%% 前五个符号发送波形和接收波形
figure;
plot((delta_t:delta_t:5*T)', x0_t(1:round(5*T/delta_t)), 'LineWidth', 1);
hold on
plot((delta_t:delta_t:5*T)', x_t(1:round(5*T/delta_t)), 'LineWidth', 1);
xlabel('时间');
ylabel('波形');
title('前5个符号波形');
legend('基带脉冲成形信号','已调信号');
box on; grid on;

y_t = x_t + noise;
figure;
plot((delta_t:delta_t:5*T)', y_t(1:round(5*T/delta_t)), 'LineWidth', 1);
xlabel('时间');
ylabel('波形');
title('接收信号波形');
box on; grid on;


% 绘制功率谱
figure;
plot(f, (fftshift(psdEstimate)));
hold on
plot(f, (S_X));
xlabel('频率 (Hz)');
ylabel('功率谱密度 ');
title('功率谱密度');
%% 

% M=4，由bit序列映射到电平符号
function mod_data = my_gray_map_real_M4(bit_data)
if mod(length(bit_data),2)==0
    % 对应的符号序列长度
    N = length(bit_data)/2;
    mod_data = zeros(1,N);
    % 逐个符号判断映射关系
    % 表1映射：10->-3A, 01->-A, 11->A, 00->3A
    for n = 1:N
        current_bits = bit_data(2*n-1 : 2*n);
        if current_bits(1)==0 && current_bits(2)==0
            mod_data(n) = 3;
        elseif current_bits(1)==0 && current_bits(2)==1
            mod_data(n) = -1;
        elseif current_bits(1)==1 && current_bits(2)==1
            mod_data(n) = 1;
        else % 10
            mod_data(n) = -3;
        end
    end
else
    error('Invalid modulation order');
end
end


%M=4，由电平符号映射到bit序列
function demod_bit_data = my_inverse_gray_map_real_M4(rx_data)
% 符号序列长度
N = length(rx_data);
demod_bit_data = zeros(1,2*N);
% 逐个符号判决，基于归一化电平 {-3, -1, 1, 3}
for n = 1:N
    current_level = rx_data(n);
    % 判决门限 -2, 0, 2
    if current_level < -2
        % -3 -> 10
        demod_bit_data(2*n-1:2*n) = [1 0];
    elseif current_level < 0
        % -1 -> 01
        demod_bit_data(2*n-1:2*n) = [0 1];
    elseif current_level < 2
        % 1 -> 11
        demod_bit_data(2*n-1:2*n) = [1 1];
    else
        % 3 -> 00
        demod_bit_data(2*n-1:2*n) = [0 0];
    end
end
end


function y_conv = conv_no_delay(y, h, delta_t)
% y: 输入信号
% h: 卷积核（匹配滤波模板）
% delta_t: 采样间隔
% y_conv: 卷积输出，长度与 y 相同，且对齐原始信号

    % 确保行向量
    y = y(:).';
    h = h(:).';
    
    Ny = length(y);
    Nh = length(h);
    
    % 初始化输出
    y_conv = zeros(1, Ny);
    
    % 直接计算卷积
    for n = 1:Ny
        % 卷积窗口
        k_min = max(1, n-Nh+1);
        k_max = n;
        h_start = Nh - (k_max - k_min);
        h_end = Nh;
        
        y_conv(n) = sum(y(k_min:k_max) .* h(h_start:h_end)) * delta_t;
    end
end