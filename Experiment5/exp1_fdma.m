% 通信与网络 实验5 FDMA
% 载波正交 BPSK
clearvars;
close all;
clc;
% rng(2023);  %固定随机种子，保证结果可复现

% 脉冲成形波形参数
T = 0.01; % 周期
V = 1/sqrt(T); % 幅度

% 用户数（>=2）
N_u = 4;

% 载波参数
delta_fc = 200; % 载波频率间隔
f_lc = 200; % 最低载波频率
f_c = f_lc + (0:N_u-1)*delta_fc; % 载波频率，满足f_c(n)*T为整数

% 噪声的单边功率谱密度
n0 = 0.5;
snr_range = -5:1:5;

% 最小时间间隔
delta_t = 1e-4;

% 每个符号周期采样次数
num_t = round(T/delta_t);

% 功率谱计算选项：1为计算，0为不计算
if_ps = 0;

N_s = 1e5; % 符号个数
N_b = N_s; % 每用户信息比特数

Es = V^2*T; % 符号能量
Eb = Es; % BPSK的比特能量

% 仿真时间
t = (delta_t:delta_t:N_s*T);

% 功率谱滑动窗参数
win_size = 100; % 滑动窗长与符号周期的比值 
win_step = 80; % 滑动窗移动步长与符号周期的比值
win_num = floor((N_s-win_size)/win_step)+1; % 滑动窗个数

% 功率谱频率分辨率
delta_f = 1/T/win_size;

% 功率谱采样频点
f = (0:win_size*num_t-1) * delta_f; % FFT所得功率谱循环周期为1/delta_t

% 成形脉冲
p_t = V * ones(num_t, 1);

% 统计ber
avg_ber_snr = zeros(1, length(snr_range));
ber_theory_snr = zeros(1, length(snr_range));

for snr_cnt = 1:length(snr_range)
    disp(['仿真Eb/n0 = ' num2str(snr_range(snr_cnt)) ' dB...']);
    if snr_cnt == length(snr_range)
        if_ps = 1; 
    end
    tic;

    % 基带脉冲成形
    x0_t = zeros(N_u, N_s*num_t);

    % 发送波形
    x_t = zeros(1, N_s*num_t);

    snr = snr_range(snr_cnt);
    n0 = Eb/10^(snr/10);

    % 随机生成符号序列
    bit_data = randi([0 1], N_u, N_b);

    % 转化为BPSK调制后的电平序列（mod_data）
    mod_data = ma_gray_map_real_M2(bit_data);

    % 脉冲成形
    for n_u = 1:N_u
        for n_s = 1:N_s
            x0_t(n_u,num_t*(n_s-1)+1:num_t*n_s) = mod_data(n_u,n_s)*p_t;
        end
    end
    % 逐用户上变频
    for n_u = 1:N_u
        x_t = x_t + x0_t(n_u, :) .* sqrt(2) .* cos(2 * pi * f_c(n_u) * t);
    end
    if if_ps == 1
        S_X = zeros(1,win_size*num_t); % 功率谱初始化
        for win_id = 1:win_num % 滑动窗编号
            % 滑动窗的时间范围
            window = (win_id-1)*win_step*num_t+1:(win_id-1)*win_step*num_t+win_size*num_t;
            % 累加各窗内波形FFT的模方
            S_X = S_X + abs(fft(x_t(window))).^2;
        end
        S_X = S_X / win_num; % 窗之间求平均
        S_X = S_X / (num_t*win_size)^2 / delta_f; % 功率谱密度函数系数
    end

    % 加性噪声
    noise_power = n0/2/delta_t;
    noise = sqrt(noise_power) * randn(size(t));

    % AWGN信道
    y_t = x_t + noise;

    % 等效基带接收波形
    y0_t = zeros(N_u, N_s*num_t);

    % 逐用户下变频
    for n_u = 1:N_u
        y0_t(n_u,:) = y_t .* sqrt(2) .* cos(2 * pi * f_c(n_u) * t);
    end

    % 初始化接收符号
    y = zeros(N_u, N_s);

    % 逐用户最佳接收
    for n_u = 1:N_u
        for n_s = 1:N_s 
            idx_start = (n_s - 1) * num_t + 1;
            idx_end = n_s * num_t;
            
            y(n_u,n_s) = sum(y0_t(n_u, idx_start:idx_end)) * delta_t;
        end
    end

    % 最佳判决结果
    bit_hat = ma_inverse_gray_map_real_M2(y);

    % 差错比特
    bit_err = xor(bit_data, bit_hat);
    % 计算各用户的误比特率
    ber = sum(bit_err,2);
    ber = ber/N_b;

    ber_theory=qfunc(sqrt(Eb/(n0/2)));

    avg_ber_snr(snr_cnt) = mean(ber);
    ber_theory_snr(snr_cnt) = ber_theory;

    toc;

    % 绘图
    fig_wave = figure('Visible', 'off');
    plot((delta_t:delta_t:4*T)', y_t(1:4*num_t), ...
            'LineWidth', 1);
    hold on;
    plot((+delta_t:delta_t:4*T)', x_t(1:4*num_t), ...
            'LineWidth', 1);
    xlabel('时间/s');
    ylabel('幅度/V');
    title(['波形图 (E_b/n_0 = ' num2str(snr) ' dB)']);
    box on;
    grid on;
    legend('接收波形信号','发送波形信号','fontsize',12);
    saveas(fig_wave, sprintf('report/images/Waveform_SNR_%gdB.png', snr));
    close(fig_wave);

    if if_ps == 1
        fig_ps = figure('Visible', 'off');
        plot(f(f<1/2/delta_t), 10*log10(S_X(f<1/2/delta_t)), 'b', 'LineWidth', 1);
        hold on;
        % FFT所得功率谱循环周期为1/delta_t，将高频部分循环位移至负频率
        plot(f(f>=1/2/delta_t)-1/delta_t, 10*log10(S_X(f>=1/2/delta_t)), 'b', 'LineWidth', 1);
        xlabel('频率/Hz');
        ylabel('功率谱密度/[dBW/Hz]');
        title('功率谱');
        xlim([-1/2/delta_t 1/2/delta_t])
        ylim(10*log10(max(S_X)*[1E-5 1.2]));
        box on;
        grid on;
        saveas(fig_ps, 'report/images/Power_Spectrum_SNR.png');
        close(fig_ps);

        % 计算总功率
        P_total = sum(S_X) * delta_f;
        disp(['总功率 ' num2str(P_total) ' W']);
        % 计算各用户功率
        for n_u = 1:N_u
            idx_range = find(abs(f - f_c(n_u)) <= delta_fc/2);
            P_user = sum(S_X(idx_range)) * delta_f;
            disp(['用户' num2str(n_u) '：中心频率 ' num2str(f_c(n_u)) ' Hz，功率 ' num2str(P_user) ' W']);
        end

    end
end

% 绘图 ber-Eb/n0
fig_ber = figure('Visible', 'off');
semilogy(snr_range, avg_ber_snr, 'b-o', 'LineWidth', 1, 'MarkerSize', 6);
hold on;
semilogy(snr_range, ber_theory_snr, 'r--x', 'LineWidth', 1, 'MarkerSize', 6);
grid on;
xlabel('E_b/n_0 (dB)');
ylabel('BER');
legend('仿真', '理论', 'Location', 'southwest');
title('用户平均 BER 曲线');
saveas(fig_ber, 'report/images/BER_Curve_Comparison.png');
close(fig_ber);


% M=2，由bit序列映射到电平符号
function mod_data = ma_gray_map_real_M2(bit_data)
% 对应的符号序列长度
[N_u, N_s] = size(bit_data);
mod_data = zeros(N_u, N_s);
% 逐用户、逐符号判断映射关系
for n_u = 1:N_u
    for n_s = 1:N_s
        current_bits = bit_data(n_u,n_s);
        if current_bits == 0
            mod_data(n_u,n_s) = -1;
        elseif current_bits == 1
            mod_data(n_u,n_s) = 1;
        end
    end
end
end

% M=2，由电平符号映射到bit序列
function demod_bit_data = ma_inverse_gray_map_real_M2(rx_data)
% 符号序列长度
[N_u, N_s] = size(rx_data);
demod_bit_data = zeros(N_u,N_s);
% 逐用户、逐符号判决
for n_u = 1:N_u
    for n_s = 1:N_s
        current_level = rx_data(n_u,n_s);
        if current_level < 0
            demod_bit_data(n_u,n_s) = 0;
        else
            demod_bit_data(n_u,n_s) = 1;
        end
    end
end
end
