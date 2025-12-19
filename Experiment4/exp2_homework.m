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
Q = $TODO$;
G = [eye(k) Q];

% 监督矩阵
H = $TODO$;

% 随机生成Mk个需要传输的比特
x_data = randi([0 1], 1, M*k);
% fprintf('x_data\n');
% disp(x_data);

% 每4个信息比特作为一个比特块，编码成（7，4）汉明码，得到比特流向量
x_code = zeros(1, M*n);
for i = 1:M
    x_code((i-1)*n+1:i*n) = mod(x_data((i-1)*k+1:i*k)*G , 2); 
end
% fprintf('x_code\n');
% disp(x_code);

% 交织 按行写入 按列读出成比特流向量
x_interleave = zeros(1,M*n);
for i = 1:number
    x_interleave((i-1)*row*column+1:i*row*column) =  interleaver(row, column, x_code((i-1)*row*column+1:i*row*column));
end
% fprintf('x_interleave\n');
% disp(x_interleave);

% 经过信道 产生长度为L的突发错误
L = 5;
y = zeros(1,M*n);
for i = 1:number
    y((i-1)*row*column+1:i*row*column) =  burst_error(x_interleave((i-1)*row*column+1:i*row*column),L);
end
% fprintf('y\n');
% disp(y);

% 解交织 按列写入 按行读出成比特流向量
y_deinterleave = zeros(1,M*n);
for i = 1:number
    y_deinterleave((i-1)*row*column+1:i*row*column) =  deinterleaver(row, column, y((i-1)*row*column+1:i*row*column));
end
% fprintf('y_deinterleave\n');
% disp(y_deinterleave);

% 每7个比特为一组，进行译码
y_decode = zeros(1,M*k);
for i = 1:M
    % 利用监督矩阵计算校正子
    syndrome = $TODO$;

    % 比较校正子和监督矩阵，找出错误位置
    if ismember(H',syndrome,'rows') == zeros(n,1)
        error_positions = 0; 
    else
        error_positions = find(ismember(H',syndrome,'rows'));
    end

    % 进行纠错
    y_decode_all = y_deinterleave((i-1)*n+1:i*n);
    if error_positions ~= 0
        y_decode_all(error_positions) = ~y_decode_all(error_positions); 
    end

    % 去除监督位
    y_decode((i-1)*k+1:i*k) = y_decode_all(1:k);
end
% fprintf('y_decode\n');
% disp(y_decode);

% 计算有交织时的误块率和误比特率
result = mod(x_data+y_decode,2);
result = transpose(reshape(result,[k,M]));
BlockErrorRate_code = sum(~ismember(result,zeros(1,k),'rows'))/M;
BitErrorRate_code = sum(result,'all')/(M*k);

function x_interleave = interleaver(row, column, x)
x = $TODO$;
x = x';
x_interleave = x(:)';
end

function y = burst_error(x, L)
noise = zeros(1,size(x,2));
error_idx = randi([1,size(x,2)-L+1]);
noise(1,error_idx:error_idx+L-1) = (rand(1,L)<0.5);
y = mod(x + noise, 2);
end

function y_deinterleave = deinterleaver(row, column, y)
y = $TODO$;
y = y';
y_deinterleave = y(:)';
end


