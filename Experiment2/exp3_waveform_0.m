% 通信与网络 实验2 基带传输
% 波形传输
clearvars;
close all;
clc;
% rng(2025);  %固定种子

% 波形参数
V = 2; % 幅度
T = 1:1:10; % 持续时间

% 相关波形的增益
alpha = 3;

% 噪声的单边功率谱密度
n0 = 8;

% 采样时间间隔
delta_t = [0.1; 0.01];

% 重复实验次数，可修改
N = 1e5;

% 符号电平均方为Ex2,电平集合为{-1,1}，等概选取
Ex2 = 1;
% 判决电平的理论噪声均方
noise_energy_theory = zeros(length(delta_t),length(T));
% 判决电平的理论信号均方
signal_energy_theory = zeros(length(delta_t),length(T));
% 判决电平的统计噪声均方
noise_energy = zeros(length(delta_t),length(T));
% 判决电平的统计信号均方
signal_energy = zeros(length(delta_t),length(T));
% 误符号率
ser = zeros(length(delta_t),length(T));
% 接收滤波器相关系数
rho = zeros(1,length(T));

% 用来画采样电平的目标T
T_rec = 10;

for cnt1 = 1:length(delta_t)
    for cnt2 = 1:length(T)
        t = (0:floor( T(cnt2) / delta_t(cnt1) )-1)'*delta_t(cnt1);
        % 成形脉冲
        p_t = V * ones(size(t));
        g_t = alpha * p_t;
        noise_energy_theory(cnt1,cnt2) = n0/2*sum(g_t.*g_t)*delta_t(cnt1);
        signal_energy_theory(cnt1,cnt2) = Ex2*(sum(p_t.*g_t)*delta_t(cnt1))^2;
        for cnt3 = 1:N
            % 随机生成一个符号电平{-1，1}
            x = $TODO$;
            % 生成一个符号的波形采样序列
            x_t = x.*p_t;
            % 加性噪声
            noise_power = $TODO$;
            noise = $TODO$;
            % AWGN信道
            y_t = x_t + noise;
            % 画采样电平
            if $TODO$ && ismember(T(cnt2),T_rec)
                figure;
                plot(t,y_t);
                box on; grid on;
                title(['T=',num2str(T_rec),', \Delta t=',num2str(delta_t(cnt1))]);
            end
            % 相关接收
            y = $TODO$;
            % 判决结果
            x_hat = $TODO$;
            % 累计误符号率
            ser(cnt1,cnt2) = ser(cnt1,cnt2) + (x~=x_hat);
            % 累计噪声能量
            noise_energy(cnt1,cnt2) = ...
                noise_energy(cnt1,cnt2) + $TODO$;
            % 累计符号能量
            signal_energy(cnt1,cnt2) = ...
                signal_energy(cnt1,cnt2) + $TODO$;
        end
    end
end

ser = ser/N;
noise_energy = noise_energy/N;
signal_energy = signal_energy/N;
% 判决电平统计信噪比
snr_dB = 10*log10(signal_energy./noise_energy);
% 判决电平理论信噪比
snr_dB_theory = 10*log10(signal_energy_theory./noise_energy_theory);

% 计算波形的理论信噪比Es = V^2*T
Es_n0_wave = Ex2*V^2*T/(n0/2);
Es_n0_wave_dB = 10*log10(Es_n0_wave);
ser_theory = $TODO$;

% 绘图

% 信噪比-T
figure; hold on;
plot(T, Es_n0_wave_dB, 'g', 'LineWidth', 3);
plot(T, snr_dB(1,:), 'bo', ...
    T, snr_dB(2,:), 'rx', ...
    T, snr_dB_theory(1,:), 'm-', ...
    T, snr_dB_theory(2,:), 'k--', ...
    'LineWidth', 1);
xlabel('波形持续时间T');
ylabel('信噪比E_s/(n_0/2)(dB)');
legend('波形的E_s/(n_0/2)', ...
    '判决电平信噪比 \Deltat=0.1, simulation', ...
    '判决电平信噪比 \Deltat=0.01, simulation', ...
    '判决电平信噪比 \Deltat=0.1, theory', ...
    '判决电平信噪比 \Deltat=0.01, theory', ...
    'Location', 'northwest');
title('信噪比-T');
box on;
grid on;

% SER-信噪比
figure;
semilogy(Es_n0_wave_dB, ser_theory, 'k', ...
    Es_n0_wave_dB, ser(1,:), 'bo', ...
    Es_n0_wave_dB, ser(2,:), 'rx', ...
    'LineWidth', 1);
xlabel('信噪比E_s/(n_0/2)(dB)');
ylabel('误符号率');
legend('Theory', ...
    '\Deltat=0.1, simulation', ...
    '\Deltat=0.01, simulation');
title('波形信道 SER-E_s/n_0');
box on;
grid on;

