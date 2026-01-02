clc; clear; close all;

rng(093);

% 码长
n = 7;
% 信息位长
k = 4;
% 传输比特块数
M = 10000;

% 交织块行数
row = 5;
% 交织块列数
column = 7;
% 交织块数量
number = (M*n)/(row*column);

% 生成矩阵
% 选取满足汉明码要求的Q矩阵
Q = [1 1 0; 0 1 1; 1 1 1; 1 0 1];
G = [eye(k) Q];

% 监督矩阵
H = [Q.' eye(n-k)];

% 随机生成Mk个需要传输的比特
x_data = randi([0 1], 1, M*k);
% fprintf('x_data\n');
% disp(x_data);

% 每4个信息比特作为一个比特块，编码成（7，4）汉明码，得到比特流向量
x_code = zeros(1, M*n);
for i = 1:M
    x_code((i-1)*n+1:i*n) = mod(x_data((i-1)*k+1:i*k)*G , 2); 
end

% 交织 按行写入 按列读出成比特流向量
x_interleave = zeros(1,M*n);
for i = 1:number
    x_interleave((i-1)*row*column+1:i*row*column) = interleaver(row, column, x_code((i-1)*row*column+1:i*row*column));
end
% 不交织 直接赋值
x_no_interleave = x_code;

% 经过信道 产生长度为L的突发错误
L_candidates = [3, 5, 10, 15, 20, 25];
fprintf("突发错误长度\t有交织误块率\t有交织误比特率\t无交织误块率\t无交织误比特率\n");
fprintf("---------------------------------------------------------------\n");
for idx1 = 1:length(L_candidates)
    L = L_candidates(idx1);
    y_interleave = zeros(1,M*n);
    y_no_interleave = zeros(1,M*n);
    for i = 1:number
        y_interleave((i-1)*row*column+1:i*row*column) =  burst_error(x_interleave((i-1)*row*column+1:i*row*column),L);
        y_no_interleave((i-1)*row*column+1:i*row*column) =  burst_error(x_no_interleave((i-1)*row*column+1:i*row*column),L);
    end

    % 解交织 按列写入 按行读出成比特流向量
    y_deinterleave = zeros(1,M*n);
    for i = 1:number
        y_deinterleave((i-1)*row*column+1:i*row*column) =  deinterleaver(row, column, y_interleave((i-1)*row*column+1:i*row*column));
    end
    % 不解交织 直接赋值
    y_no_deinterleave = y_no_interleave;

    % 每7个比特为一组，进行译码
    y_decode_interleave = zeros(1,M*k);
    y_decode_no_interleave = zeros(1,M*k);
    for i = 1:M
        % 利用监督矩阵计算校正子
        syndrome_interleave = mod(y_deinterleave((i-1)*n+1:i*n) * H', 2);
        syndrome_no_interleave = mod(y_no_deinterleave((i-1)*n+1:i*n) * H', 2);
        % 比较校正子和监督矩阵，找出错误位置
        if ismember(H',syndrome_interleave,'rows') == zeros(n,1)
            error_positions_interleave = 0; 
        else
            error_positions_interleave = find(ismember(H',syndrome_interleave,'rows'));
        end
        if ismember(H',syndrome_no_interleave,'rows') == zeros(n,1)
            error_positions_no_interleave = 0; 
        else
            error_positions_no_interleave = find(ismember(H',syndrome_no_interleave,'rows'));
        end
        % 进行纠错
        y_decode_all = y_deinterleave((i-1)*n+1:i*n);
        y_decode_no_all = y_no_deinterleave((i-1)*n+1:i*n);
        if error_positions_interleave ~= 0
            y_decode_all(error_positions_interleave) = ~y_decode_all(error_positions_interleave); 
        end
        if error_positions_no_interleave ~= 0
            y_decode_no_all(error_positions_no_interleave) = ~y_decode_no_all(error_positions_no_interleave); 
        end
        % 去除监督位
        y_decode_interleave((i-1)*k+1:i*k) = y_decode_all(1:k);
        y_decode_no_interleave((i-1)*k+1:i*k) = y_decode_no_all(1:k);
    end

    % 计算有交织时的误块率和误比特率
    result_interleave = mod(x_data+y_decode_interleave,2);
    result_interleave = transpose(reshape(result_interleave,[k,M]));
    BlockErrorRate_code_interleave = sum(~ismember(result_interleave,zeros(1,k),'rows'))/M;
    BitErrorRate_code_interleave = sum(result_interleave,'all')/(M*k);
    % 计算无交织时的误块率和误比特率
    result_no_interleave = mod(x_data+y_decode_no_interleave,2);
    result_no_interleave = transpose(reshape(result_no_interleave,[k,M]));
    BlockErrorRate_no_interleave = sum(~ismember(result_no_interleave,zeros(1,k),'rows'))/M;
    BitErrorRate_no_interleave = sum(result_no_interleave,'all')/(M*k);

    fprintf("%d\t%.5f\t%.5f\t%.5f\t%.5f\n", L, BlockErrorRate_code_interleave, BitErrorRate_code_interleave, BlockErrorRate_no_interleave, BitErrorRate_no_interleave);
end
