clc
clear
% 图拓扑加载
load('topology_general.mat');
% 获取节点数量
N = max([s; t], [], 'all');
% 选择
viewNode = 1; % TODO 1<=viewNode<=N,选中某一路由节点，观察该路由节点的最小生成树随链路状态传播而变化的情况
% 路由节点编号 "a"->1 "b"->2, ... ; 共N个编号
node = 1:N;
% 获取连边数量
E = length(s);
% 邻接图，值为代价，inf代表不连通
G = inf(N, N);
for i = 1:E
    G(s(i), t(i)) = w(i);
    G(t(i), s(i)) = w(i);
end

% 链路状态法（Dijkstra 算法）
seenGs = inf(N, N, N);  % 每个节点可见的链路状态（源节点×目标节点×每个节点）
for i = 1:N
    seenGs(i, :, i) = G(i, :);  % 初始化，每个节点知道自己（源）和自己邻居（目标）的链路状态
end

% 假设链路状态在开始后一起进行泛洪，每次循环传输一跳，每次循环后，各个节点根据已知的链路状态计算路由表
% 路由表体现为最小生成树
% 某一节点s发的LSA通告为G(s, :)中非inf的元素

% 从s发出的链路状态包的大小（不含到自身的0代价"链路"）
LSASize = zeros(1,N);
for i = 1:N
    LSASize(i) = sum(G(i, :) ~= inf) - 1;  
end

% 定义两个路由器中的缓冲区，记录收到的LSA通告
% bufferLast(i,j)={k1,k2}，表示上一轮传输中，路由节点i收到了节点j转发的源自路由节点k1,k2发出的LSA通告; 否则, buffer_L(i,j)=0
bufferLast = cell(N,N);
% buffer(i,j)={k1,k2}，表示本轮传输中，路由节点i将收到了节点j转发的源自路由节点k1,k2发出的LSA通告;否则, buffer(i,j)=0
buffer = cell(N,N);
% 接收到的LSA通告，received(i,j)=true表示路由节点i收到过源自路由节点k发出的LSA通告
received = false(N,N);

% 开始模拟泛洪
transferCount = 0;  % 总传输次数
transferSize = 0;  % 总传输包大小（每条路由记录认为是大小1）
% 初始泛洪
currentCount = 0;
currentSize = 0;
for i = 1:N
    neighbors = node(G(i, :) ~= inf);  % 节点s的邻居节点
    % 向邻居节点泛洪LSA
    LSASource = i;
    for n = neighbors
        % 向cell中添加一个值
        buffer{n,i}(end+1) = LSASource;  
        % 当前轮次泛洪数据统计
        currentCount = currentCount + 1;
        currentSize = currentSize + LSASize(LSASource);
    end
end
transferCount = transferCount + currentCount;
transferSize = transferSize + currentSize;
% 后续的泛洪过程
numLoop = 100; % 最大迭代轮次
for loop = 1:numLoop
    % 刷新buffer
    bufferLast = buffer;
    buffer = cell(N,N);

    % 统计
    currentCount = 0;
    currentSize = 0;

    for i = 1:N % 路由节点s处理bufferLast中收到的LSA通告
        for j = 1:N  % 处理每个相邻路由器发来的LSA
            for LSASource = bufferLast{i, j}  % 逐个LSA处理
                if received(i, LSASource) % 该LSA通告已经收到过
                    continue;
                end
                % 节点i更新LSA数据库
                seenGs(LSASource, :, i) = G(LSASource,:);
                seenGs(:, LSASource, i) = G(:,LSASource);
                received(i, LSASource) = true;
                % 泛洪
                neighbors = node(G(i, :) ~= inf);
                for n = neighbors
                    if n == j  % 不向 向自己转发LSA的路由器 泛洪
                        continue
                    end
                    buffer{n,i}(end+1) = LSASource;  % 向cell中添加一个值
                    currentCount = currentCount + 1;
                    currentSize = currentSize + LSASize(LSASource);
                end
            end
        end
    end
    fprintf('完成第%d次泛洪迭代...\n', loop);
 
    %统计数据更新
    transferCount = transferCount + currentCount;
    transferSize = transferSize + currentSize;

    % 每个节点采用dijkstra算法计算最小生成树
    for i = 1:N
        MST = dijkstra(seenGs(:, :, i), i);
        if i == viewNode
            MST
            figure;
            mstG = graph(MST, nodeNames);
            plot(mstG, 'XData', nodePositions(:, 1), 'YData', nodePositions(:, 2), 'EdgeLabel', mstG.Edges.Weight);
            title(sprintf('节点%d在收到第%d轮通过泛洪到达的链路状态后的最小生成树', i, loop));
            axis equal;
            hold on;
        end
    end

    % 结束的判断
    if currentCount == 0 %各节点LSA数据库已经收敛
        fprintf('路由收敛，总传输次数 %d 次，总传输包大小 %d。\n', transferCount, transferSize);
        break;
    end
    if loop == numLoop
        fprintf('达到循环次数上限，路由信息计算结束。\n');
    end
end

function MST = dijkstra(G, s)
    % 图中的节点数
    N = size(G, 1);
    % 初始化数组
    dist = inf(1, N);   % 从起始节点到每个节点的距离
    visited = false(1, N);  % 跟踪已访问的节点
    parent = zeros(1, N);   % 最短路径中的父节点
    % 起始节点到自身的距离为0
    dist(s) = 0;
    for i = 1:N-1
        % 找到距离最小的节点
        u = findMinDistance(dist, visited);
        if u == -1  % 没有连通的未访问节点
            break
        end
        % 将节点标记为已访问
        visited(u) = true;
        % 更新相邻节点的距离
        for v = 1:N
            if isinf(G(u, v))  % 不相邻
                continue
            end
            % 对于未被访问的节点，更新到源点的距离
            if ~visited(v) && dist(u) + G(u, v) < dist(v)
                dist(v) = dist(u) + G(u, v);
                parent(v) = u;
            end
        end
    end

    % 构建最小生成树
    MST = zeros(N, N);
    for v = 1:N
        if v ~= s && parent(v) ~= 0
            MST(v, parent(v)) = G(v, parent(v));
            MST(parent(v), v) = G(parent(v), v);
        end
    end
end

% 找到距离最小的节点
function u = findMinDistance(dist, visited)
    minDist = inf;
    u = -1;
    for i = 1:length(dist)
        if ~visited(i) && dist(i) < minDist
            minDist = dist(i);
            u = i;
        end
    end
end