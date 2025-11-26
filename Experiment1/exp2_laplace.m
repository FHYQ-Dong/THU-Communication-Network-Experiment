clc; clear; close all;
rng(2025);  % 固定随机种子，保证结果可复现


%% -----  参数  -----
V_max = 10;
sigma_x = sqrt(2);  % 拉普拉斯分布的系数
N = 1000000;  % 拉普拉斯分布下采样点数
A = randlap(N, sigma_x);  % 拉普拉斯分布下的采样
mu = 255; % μ取值
% 画出幅度分布
bins = -V_max:2*V_max/100:V_max; % 选取合理的直方图格点数
histogram(A, bins,'Normalization','probability');
title("点列A的幅度分布");


%% -----  A均匀量化，注意对信号截断到[-10,10]区间内  -----
A_hat = zeros(size(A));
L = 2^8; % 量化区间
interval = 2 * V_max / L; % 量化间隔
x = -V_max:interval:V_max; % 分层电平
y = (x(1:end-1) + x(2:end)) / 2; % 重建电平
for i = 1:N  % 循环点列
    for j = 1:L  % 循环区间
        % 完成均匀量化，得到量化后的采样A_hat
        if A(i) >= x(j) && A(i) < x(j+1)
            A_hat(i) = y(j);
            break;
        elseif A(i) <= x(1) % 小于最小值
            A_hat(i) = y(1);
            break;
        elseif A(i) >= x(end) % 大于最大值
            A_hat(i) = y(end);
            break;
        end
    end
end
ex = A - A_hat; % 均匀量化的量化误差e(x)
power_exp = sum(A.^2) / N;  % 信号功率
noise_exp = sum(ex.^2) / N; % 噪声方差
snr_exp = power_exp/noise_exp;
snrdb_exp = 10*log10(snr_exp);

fprintf("均匀量化下，量化步长=%f\n",interval);
figure;
nbins = 100; % 选取合理的直方图格点数
histogram(ex,nbins,'Normalization','probability');
title(sprintf("均匀量化的量化误差分布（区间数：%d）", nbins));
xlabel("e"); ylabel("p(e)");
figure;
histogram(ex, -interval:interval/100:interval,'Normalization','probability');
xlabel("e"); ylabel("p(e)");
title("均匀量化的量化误差分布（步长视图，区间数：100）");
fprintf("均匀量化下，量化噪声方差%f,点列A信号功率%f,量化信噪比%fdB\n",noise_exp,power_exp,snrdb_exp);


%% μ律压缩
B = mu_compress(A, V_max, mu);
figure;
nbins_B = 100; % 选取合理的直方图格点数
histogram(B,nbins_B,'Normalization','probability');
title("μ律压缩后的样本分布")
xlabel("x")
ylabel("p(x)")


%% BC均匀量化，注意对信号截断到[-10,10]区间内
C = zeros(size(B));
L = 2^8; % 量化区间
interval = 2 * V_max / L; % 量化间隔
x = -V_max:interval:V_max; % 分层电平
y = (x(1:end-1) + x(2:end)) / 2; % 重建电平
for i = 1:N  % 循环点列
    for j = 1:L  % 循环区间
        % 完成均匀量化，得到量化后的采样C
        if B(i) >= x(j) && B(i) < x(j+1)
            C(i) = y(j);
            break;
        elseif B(i) <= x(1) % 小于最小值
            C(i) = y(1);
            break;
        elseif B(i) >= x(end) % 大于最大值
            C(i) = y(end);
            break;
        end
    end
end
ex1 = B - C; % BC间的量化误差e(x)
powerBC_exp = sum(B.^2) / N;  % 信号功率
noiseBC_exp = sum(ex1.^2) / N; % 噪声方差
snrBC_exp = powerBC_exp/noiseBC_exp;
snrBCdb_exp = 10*log10(snrBC_exp);
figure;
nbins_ex1 = 100; % 选取合理的直方图格点数
histogram(ex1,nbins_ex1,'Normalization','probability');
xlabel("e"); ylabel("p(e)");
title(sprintf("BC间的量化误差分布（区间数：%d）", nbins_ex1));
figure;
histogram(ex1, -interval:interval/100:interval,'Normalization','probability');
xlabel("e")
ylabel("p(e)")
title("BC间的量化误差分布（步长视图，区间数：100）");
fprintf("压缩后均匀量化下，量化步长=%f\n",interval);
fprintf("BC间,量化噪声方差%f,点列B信号功率%f,量化信噪比%fdB\n",noiseBC_exp,powerBC_exp,snrBCdb_exp);


%% μ律扩张
D = mu_expand(C, V_max, mu);
ex2 = A - D; % AD间的量化误差e(x)
powerAD_exp = sum(A.^2) / N;  % 信号功率
noiseAD_exp = sum(ex2.^2) / N; % 噪声方差
snrAD_exp = powerAD_exp/noiseAD_exp;
snrADdb_exp = 10*log10(snrAD_exp);
figure;
nbins_ex2 = 100; % 选取合理的直方图格点数
histogram(ex2,nbins_ex2,'Normalization','probability');
title(sprintf("AD间的量化误差分布（区间数：%d）", nbins_ex2));
xlabel("e")
ylabel("p(e)")
fprintf("AD间,量化噪声方差%f,点列A信号功率%f,量化信噪比%fdB\n",noiseAD_exp,powerAD_exp,snrADdb_exp);
figure;
histogram(ex2, -interval:interval/100:interval,'Normalization','probability');
xlabel("e")
ylabel("p(e)")
title("AD间的量化误差分布（步长视图，区间数：100）");


%% 采用逆累积分布函数生成零均值拉普拉斯采样
function x = randlap(siz, sigma_x)
    x = (log(rand(siz,1)).*(2*floor(rand(siz,1)*2)-1))*sigma_x/sqrt(2);
end

