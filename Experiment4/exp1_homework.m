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
H = $TODO$;

% 随机生成M个需要传输的比特块
x_data = randi([0 1], M, k);

% 将每个比特块编码成（7，4）汉明码
x_code = mod(x_data*G , 2); %进行编码


% 经过误符号率为p的BSC信道
p = 0.001;
noise = rand(M, n) < p;
y = mod(x_code + noise, 2);
 
% 利用监督矩阵计算校正子
syndrome = $TODO$;
 
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
result_code = $TODO$;
BlockErrorRate_code = $TODO$;
BitErrorRate_code = $TODO$;
