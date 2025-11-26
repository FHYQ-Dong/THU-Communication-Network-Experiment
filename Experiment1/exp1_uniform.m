clc; clear; close all;
rng(2025);  % 固定随机种子，保证结果可复现


%% -----  参数  -----
V_max = 1;  % 量化范围[-V_max, V_max]
Q_bit = [1,2,3];  % 量化比特数
V = 1;  % 均匀分布的范围[-V, V]
N = 1000000;  % 均匀分布下采样点数
s = -V + 2*V*rand(1, N);  % 均匀分布下的采样
s_hat = zeros(size(s));  % 量化后的重建采样
% 画出幅度分布
bins = 100; % 选取合理的直方图格点数
histogram(s, bins,'Normalization','probability');
title(sprintf("点列A的幅度分布（区间数：%d）", bins));


%% -----  实验  -----
% 1bit, 2bit, 3bit均匀量化实验
for B = Q_bit
    L = 2 ^ B; % 量化区间
    interval = 2 * V / L; % 量化间隔
    x = linspace(-V, V, L+1); % 分层电平[x_1,x_2,..,x_{L+1}]
    y = (x(1:end-1) + x(2:end)) / 2; % 重建电平[y_1,y_2,...,y_L]

    % 理论计算
    noise_th = interval ^ 2 / 12; % 量化噪声方差
    power_th = V ^ 2 / 3; % 信号功率
    snr_th = power_th / noise_th;
    snrdb_th = 10 * log10(snr_th);

    % 实验计算
    for i = 1:N  % 循环点列
        for j = 1:L  % 循环区间
            if s(i) >= x(j) && s(i) < x(j+1)
                s_hat(i) = y(j);
                break;
            end
        end
    end
    ex = s - s_hat;  % 量化误差e(x)
    power_exp = sum(s.^2) / N;  % 信号功率
    noise_exp = sum(ex.^2) / N; % 量化噪声方差
    snr_exp = power_exp/noise_exp;  %量化信噪比
    snrdb_exp = 10*log10(snr_exp);

    % 输出
    figure;
    exbins = -interval:interval/100:interval;
    histogram(ex, exbins, 'Normalization', 'probability'); % 计算e(x)的分布
    title(num2str(B) + sprintf("bit下的量化误差分布（区间数：%d）", 100))
    xlabel("e")
    ylabel("p(e)")
    fprintf("%dbit量化时,噪声方差理论值/实际值=%f/%f,信号功率理论值/实际值=%f/%f,量化信噪比理论值/实际值=%f/%fdB\n",B,noise_th,noise_exp,power_th,power_exp,snrdb_th,snrdb_exp);
end
