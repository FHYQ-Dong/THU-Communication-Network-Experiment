clc; clear; close all;

rng(09393);

% 码长
n = 7;
% 信息位长
k = 4;
% 传输比特块数
M = 10000;

% 生成矩阵
Q = [1 1 1; 1 1 0; 1 0 1; 0 1 1];
G = [eye(k) Q];
% 监督矩阵
H = [Q.' eye(n-k)];
% 随机生成M个需要传输的比特块
x_data = randi([0 1], M, k);
% 将每个比特块编码成（7，4）汉明码
x_code = mod(x_data*G , 2); %进行编码

% 经过误符号率为p的BSC信道
p_candidates = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2];
fprintf("信道误符号率\t无信道编码误块率\t无信道编码误比特率\t有信道编码误块率\t有信道编码误比特率\n");
fprintf("---------------------------------------------------------------\n");
for idx1 = 1:length(p_candidates)
    p = p_candidates(idx1);
    noise = rand(M, n) < p;
    y = mod(x_code + noise, 2);
    
    % 利用监督矩阵计算校正子
    syndrome = mod(y * H.', 2);
    
    % 比较校正子和监督矩阵，找出错误位置
    error_positions = zeros(M,1);
    for i = 1:M
        if ismember(H',syndrome(i,:),'rows') == zeros(n,1)
            error_positions(i,1) = 0; 
        else
            error_positions(i,1) = find(ismember(H',syndrome(i,:),'rows'));
        end
    end

    % 进行纠错
    y_decode = y;
    for i = 1:M
        if error_positions(i,1) ~= 0
            y_decode(i, error_positions(i,1)) = ~y_decode(i, error_positions(i,1)); 
        end
    end
    % 去掉监督位
    y_decode = y_decode(:,1:k);

    % 计算无信道编码时的误块率和误比特率
    result_uncode = mod(x_data+y(:,1:k),2);
    BlockErrorRate_uncode = sum(~ismember(result_uncode,zeros(1,k),'rows'))/M;
    BitErrorRate_uncode = sum(result_uncode,'all')/(M*k);

    % 计算有信道编码时的误块率和误比特率
    result_code = mod(x_data + y_decode, 2);
    BlockErrorRate_code = sum(~ismember(result_code, zeros(1, k), 'rows')) / M;
    BitErrorRate_code = sum(result_code, 'all') / (M * k);
    fprintf("%.3f\t%.5f\t%.5f\t%.5f\t%.5f\n", p, BlockErrorRate_uncode, BitErrorRate_uncode, BlockErrorRate_code, BitErrorRate_code);
end

%% 画图像
% 经过误符号率为p的BSC信道
p_log_candidates = linspace(log10(0.001), log10(0.2), 20);
p_candidates = 10.^p_log_candidates;
BitErrorRate_uncode = zeros(1, length(p_candidates));
BlockErrorRate_uncode = zeros(1, length(p_candidates));
BlockErrorRate_code = zeros(1, length(p_candidates));
BitErrorRate_code = zeros(1, length(p_candidates));

for idx1 = 1:length(p_candidates)
    p = p_candidates(idx1);
    noise = rand(M, n) < p;
    y = mod(x_code + noise, 2);
    
    % 利用监督矩阵计算校正子
    syndrome = mod(y * H.', 2);
    
    % 比较校正子和监督矩阵，找出错误位置
    error_positions = zeros(M,1);
    for i = 1:M
        if ismember(H',syndrome(i,:),'rows') == zeros(n,1)
            error_positions(i,1) = 0; 
        else
            error_positions(i,1) = find(ismember(H',syndrome(i,:),'rows'));
        end
    end

    % 进行纠错
    y_decode = y;
    for i = 1:M
        if error_positions(i,1) ~= 0
            y_decode(i, error_positions(i,1)) = ~y_decode(i, error_positions(i,1)); 
        end
    end
    % 去掉监督位
    y_decode = y_decode(:,1:k);

    % 计算无信道编码时的误块率和误比特率
    result_uncode = mod(x_data+y(:,1:k),2);
    BlockErrorRate_uncode(idx1) = sum(~ismember(result_uncode,zeros(1,k),'rows'))/M;
    BitErrorRate_uncode(idx1) = sum(result_uncode,'all')/(M*k);

    % 计算有信道编码时的误块率和误比特率
    result_code = mod(x_data + y_decode, 2);
    BlockErrorRate_code(idx1) = sum(~ismember(result_code, zeros(1, k), 'rows')) / M;
    BitErrorRate_code(idx1) = sum(result_code, 'all') / (M * k);
end

figure; hold on;
semilogy(p_candidates, BitErrorRate_uncode, '-o', 'LineWidth', 1, 'Color', 'r');
semilogy(p_candidates, BlockErrorRate_uncode, '-+', 'LineWidth', 1, 'Color', 'r');
semilogy(p_candidates, BitErrorRate_code, '-o', 'LineWidth', 1, 'Color', 'b');
semilogy(p_candidates, BlockErrorRate_code, '-+', 'LineWidth', 1, 'Color', 'b');
legend('无信道编码误比特率', '无信道编码误块率', '有信道编码误比特率', '有信道编码误块率', 'Location', 'southeast');
xlabel('信道误符号率');
ylabel('误码率');
title('误码率与信道误符号率的关系');
set(gca, 'XScale', 'log', 'YScale', 'log');
grid on; box on;
