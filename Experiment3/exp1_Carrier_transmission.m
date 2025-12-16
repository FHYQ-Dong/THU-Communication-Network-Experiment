% 通信与网络 实验3 载波传输
%BPSK
clearvars;
close all;
clc;

N = 1e5;          % 符号个数
Es = 1;           % 每符号能量
Eb = Es;          % 对BPSK, Eb=Es
T = 0.01;         % 符号周期
A = sqrt(Es/T);   % 幅度
f_c = 500;        % 载波频率
delta_t = 1e-4;   % 采样间隔
num = T/delta_t;  % 每个符号周期采样点数
t = (delta_t:delta_t:N*T);  % 仿真时间轴


% 成形脉冲
p_t = A * ones(num,1);

%初始化接收符号、判决符号
y=zeros(1,N);
x_hat=zeros(1,N);



% 随机比特序列与调制
bit_data = randi([0 1], 1, N);
mod_data = my_gray_map_real_M2(bit_data);

% 基带信号成形
x0_t = zeros(1,N*num);
for cnt = 1:N
    x0_t(num*(cnt-1)+1:num*cnt) = mod_data(cnt)*p_t;
end



%脉冲成形
for cnt = 1:N   
    x0_t(num*(cnt-1)+1:num*cnt) = mod_data(cnt)*p_t;
end
%BPSK调制
x_t = x0_t.*(sqrt(2)*cos(2*pi*f_c*t));


EbN0_dB = 2:2:8;                  % 不同信噪比(dB)
EbN0 = 10.^(EbN0_dB/10);           % 线性比值
EsN0 = EbN0; 
EsN0_dB = 10.*log10(EsN0);
N0_list = Eb ./ EbN0 ;          % 单边噪声功率谱密度



% ====================== 匹配滤波器定义 ======================
% 方法1：带通匹配滤波器
h1_t  = sqrt(2/T) * cos(2*pi*f_c*(0:delta_t:T-delta_t));  
% 方法2：基带匹配滤波器
h2_t = (1/sqrt(T)) * ones(1,num);  

%最佳采样时刻
t_sample = T:T:N*T; 

idx_sample = round(t_sample/delta_t);


%% 不加噪声观察匹配滤波

y_t = x_t;

% % 匹配滤波输出

%方法1
mf_out1 = conv_no_delay(y_t, fliplr(h1_t), delta_t);  

%方法2
y_base = y_t .* (sqrt(2)*cos(2*pi*f_c*t));
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

offset_list = [0 1 2]; % 采样点偏离 0，1，2
for idx = 1:length(N0_list)
    n0 = N0_list(idx);

    % AWGN噪声
    noise_power = n0/2/delta_t;
    noise = sqrt(noise_power) * randn(size(t));
    
    % 接收信号
    y_t = x_t + noise;
    
    % 匹配滤波
    mf1_out = conv_no_delay(y_t, fliplr(h1_t), delta_t);
    % 方法2先下变频
    y_base = y_t .* (sqrt(2)*cos(2*pi*f_c*t));
    mf2_out = conv_no_delay(y_base, fliplr(h2_t), delta_t);
    
    % 遍历采样偏差
     for j = 1:length(offset_list)
         offset_time = offset_list(j);              
    
        idx_s = idx_sample + offset_time;    
               
        new_ind = idx_s(1:length(idx_s)-1);  %避免最后一个符号采样偏移溢出，统计1~N-1个符号
        % 判决
        x_hat1 = my_inverse_gray_map_real_M2(mf1_out(new_ind));
        x_hat2 = my_inverse_gray_map_real_M2(mf2_out(new_ind));
    
        % 计算误比特
        bit_err1 = xor(bit_data(1:N-1), x_hat1);
        bit_err2 = xor(bit_data(1:N-1), x_hat2);
    
        BER_sim1(idx,j) = sum(bit_err1)/length(bit_err1);
        BER_sim2(idx,j) = sum(bit_err2)/length(bit_err2);
    end

end

% BPSK理论BER
BER_theory = qfunc(sqrt(2*EbN0));  
SER_theory = BER_theory;            % 对BPSK, SER=BER
SER_sim1 = BER_sim1;
SER_sim2 = BER_sim2;

%% ====================== 绘制BER/SER曲线 ======================

colors = lines(length(offset_list)); % 生成直观颜色
markers = {'o','s','d'};     % 不同偏差标记
figure; hold on;

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
semilogy(EbN0_dB, BER_theory, '-', 'Color', 'r',  'LineWidth',1.5);
legend_entries{end+1} = '理论值';
set(gca, 'YScale', 'log');
legend(legend_entries,'Location','best');
xlabel('E_b/N_0 (dB)');
ylabel('误码率 (BER)');
title('BPSK误比特率性能曲线');
grid on; box on;




figure;hold on;
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

semilogy(EsN0_dB, SER_theory, '-', 'Color', 'r',  'LineWidth',1.5);
legend_entries{end+1} = '理论值';

legend(legend_entries,'Location','best');
xlabel('E_s/N_0 (dB)');
ylabel('误比特率 (BER)');
title('BPSK误码率性能曲线');
set(gca, 'YScale', 'log');
grid on; box on;


%% 随采样偏差变化的误码率
offset_list2 = -3:3; 
BER_sim1_2 = zeros(size(offset_list2));
BER_sim2_2 = zeros(size(offset_list2));
n0 = N0_list(4);

% AWGN噪声
noise_power = n0/2/delta_t;
noise = sqrt(noise_power) * randn(size(t));

% 接收信号
y_t = x_t + noise;

% 匹配滤波
mf1_out = conv_no_delay(y_t, fliplr(h1_t), delta_t);
y_base = y_t .* (sqrt(2)*cos(2*pi*f_c*t));
mf2_out = conv_no_delay(y_base, fliplr(h2_t), delta_t);

% 遍历采样偏差
for j = 1:length(offset_list2)
    offset_time = offset_list2(j);
    idx_s = idx_sample + offset_time;
    % 保护索引不越界
    valid_mask = (idx_s >= 1) & (idx_s <= length(mf1_out));
    idx_s_valid = idx_s(valid_mask);
    
    % 这里简单处理，取前N-1个点
    new_ind = idx_s(1:N-1);
    
    x_hat1 = my_inverse_gray_map_real_M2(mf1_out(new_ind));
    x_hat2 = my_inverse_gray_map_real_M2(mf2_out(new_ind));
    
    bit_err1 = xor(bit_data(1:N-1), x_hat1);
    bit_err2 = xor(bit_data(1:N-1), x_hat2);
    
    BER_sim1_2(j) = sum(bit_err1)/length(bit_err1);
    BER_sim2_2(j) = sum(bit_err2)/length(bit_err2);
end

figure; hold on;
semilogy(offset_list2, BER_sim1_2,  'LineWidth',1.5);
semilogy(offset_list2, BER_sim2_2,  'LineWidth',1.5);
legend('方法1','方法2');
xlabel('采样偏差');
ylabel('误码率');
title('误码率随采样偏差变化');
set(gca, 'YScale', 'log');
grid on; box on;


%% 发送信号功率谱
% 设置滑动窗口参数
window_length = 1024; % 窗口长度
overlap = 512; % 重叠长度
nfft = 1024; % FFT长度

% 初始化PSD估计
psdEstimate = zeros(nfft, 1);  

% 滑动窗口处理
for start = 1:(window_length - overlap):(length(x_t) - window_length)
    % 提取窗口
    windowedSegment = x_t(start:start+window_length-1);
    
    % 窗函数处理（这里使用汉明窗）
    windowedSegment = windowedSegment .* hamming(window_length)';
    
    % 计算傅里叶变换
    Y = fft(windowedSegment, nfft);
    
    % 双边 = 单边 / 2
    Y = Y / 2;
    
    % 计算功率谱密度
    P2 = abs(Y/nfft).^2 ;
    
    % 将功率谱加入估计值（包括负频率）
    psdEstimate = psdEstimate + P2';
end

% 完成平均
psdEstimate = psdEstimate / ((length(x_t) - window_length) / (window_length - overlap));

% 计算频率轴
f = (-nfft/2:nfft/2-1) / (nfft * delta_t);  
%功率谱理论值
S_X = (A^2*T/2) * ((sinc((f-f_c)*T)).^2 + (sinc((f+f_c)*T)).^2);


%% 前十个符号发送波形和接收波形
figure;
plot((delta_t:delta_t:10*T)', x0_t(1:round(10*T/delta_t)), 'LineWidth', 1);
hold on
plot((delta_t:delta_t:10*T)', x_t(1:round(10*T/delta_t)), 'LineWidth', 1);
xlabel('时间');
ylabel('波形');
title('前10个符号波形');
legend('基带脉冲成形信号','已调信号');
box on; grid on;

y_t = x_t + noise;
figure;
plot((delta_t:delta_t:10*T)', y_t(1:round(10*T/delta_t)), 'LineWidth', 1);
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
title('Averaged Power Spectrum');

%% 

% M=2，由bit序列映射到电平符号
function mod_data = my_gray_map_real_M2(bit_data)
% 对应的符号序列长度
N = length(bit_data);
mod_data = zeros(1,N);
% 逐个符号判断映射关系
for n = 1:N
    current_bits = bit_data(n);
    if current_bits == 0
        mod_data(n) = -1; % 0 映射为 -1
    else
        mod_data(n) = 1;  % 1 映射为 1
    end
end
end

% M=2，由电平符号映射到bit序列
function demod_bit_data = my_inverse_gray_map_real_M2(rx_data)
% 符号序列长度
N = length(rx_data);
demod_bit_data = zeros(1,N);
% 逐个符号判决
for n = 1:N
    current_level = rx_data(n);
    if current_level < 0
        demod_bit_data(n) = 0;
    else
        demod_bit_data(n) = 1;
    end
end
end
%时间对齐的卷积函数
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