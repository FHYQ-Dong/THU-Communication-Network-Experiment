clc
clear
% 图拓扑加载
load('topology_general.mat');
% 获取节点数量
N = max([s; t], [], 'all');
% 获取连边数量
E = length(s);
fprintf("拓扑节点数: %d，拓扑边数: %d\n", N, E);
% 邻接图，值为代价，inf代表不连通
G = inf(N, N);
for i = 1:E
    G(s(i), t(i)) = w(i);
    G(t(i), s(i)) = w(i);
end
% 初始化距离矢量矩阵和下一跳矩阵（行表示源节点，列表示目标节点）
% distMatrix(n, :)代表节点n所维护的距离矢量
distMatrix = inf(N, N);
% nextHopMatrix(n, :)代表节点n所维护的下一跳路由器表
nextHopMatrix = zeros(N, N);
for n = 1:N
    distMatrix(n, n) = 0;
end
for i = 1:E
    distMatrix(s(i), t(i)) = w(i);
    distMatrix(t(i), s(i)) = w(i);
    nextHopMatrix(s(i), t(i)) = t(i);
    nextHopMatrix(t(i), s(i)) = s(i);
end

% 迭代更新距离矩阵和下一跳矩阵（Jacobi式严格同步）
[distMatrix, nextHopMatrix] = updateRouteJacobi(G, distMatrix, nextHopMatrix, 100);

% 距离向量法（严格同步 Bellman-Ford 算法）
% 所有节点用第 k 轮的旧值计算，统一生成第 k+1 轮新值
function [distMatrix, nextHopMatrix] = updateRouteJacobi(G, distMatrix, nextHopMatrix, loops)
    % 获取节点数量
    numNodes = size(G, 1);
    loopCount = 0;
    % 输出路由表
    fprintf('--------------- LOOP START --------------\n')
    disp('距离矢量矩阵:');
    disp(distMatrix);
    disp('下一跳矩阵:')
    disp(nextHopMatrix);

    transferCount = 0;  % 总传输次数
    transferSize = 0;  % 总传输包大小（每条路由记录认为是大小1）
    
    tmpdistMatrix = distMatrix;
    tmpnextHopMatrix = nextHopMatrix;
    for numLoop = 1:loops
        updated = false; % 是否收敛
        updateID = zeros(1, numNodes); % 节点在一次循环后是否更新路由表

        % 阶段1：所有节点广播（基于第 k 轮的旧表）
        for n = 1:numNodes % 遍历所有节点（进行广播）
            % 从n发出的距离矢量（使用旧表！）
            distVector = distMatrix(n, :);
            
            % 距离矢量表大小（inf的视为不存在）
            dvSize = sum(~isinf(distVector), "all") - 1;
            for m = 1:numNodes % 遍历所有相邻路由器（接收到路由信息）
                if isinf(G(m, n)) % 不相邻
                    continue
                end

                % 第一次统计：传输次数+=1，传输包大小+=距离矢量个数
                if numLoop == 1
                    transferCount = transferCount + 1;
                    transferSize = transferSize + dvSize;
                end

                % 通过邻居的路由（n的路由信息 + m->n的代价）
                altDist = G(m, n) + distVector;
                
                % 路由更新策略：
                % 1. 如果节点m去住某些目的节点的当前最佳路由下一跳为节点n，则更新这些目的节点的距离信息为altDist中对应位置的值
                nextHopIsNAndNeedUpdate = ((nextHopMatrix(m, :) == n).*(altDist ~= distMatrix(m, :)))==1;
                tmpdistMatrix(m, nextHopIsNAndNeedUpdate) = altDist(nextHopIsNAndNeedUpdate);
                
                % 2. 如果节点m去住某些目的节点的当前最佳路由下一跳不为节点n，如果改为经过n可以获得更短距离，则将这些目的节点的最佳路由下一跳改为n，且更新相应目的节点的最短距离
                betterRoute = ((~nextHopIsNAndNeedUpdate).*(altDist < distMatrix(m, :)))==1;
                tmpdistMatrix(m, betterRoute) = altDist(betterRoute);
                tmpnextHopMatrix(m, betterRoute) = n;
                
                if sum(nextHopIsNAndNeedUpdate, 'all') > 0 || sum(betterRoute, 'all') > 0
                    updated = true;
                    updateID(n) = 1;
                end
            end
        end
        
        % 若节点的路由表更新，应该广播给邻居节点
        for n = 1:numNodes
            if updateID(n) == 1
                distMatrix = tmpdistMatrix;
                nextHopMatrix = tmpnextHopMatrix;
                dvSize = sum(~isinf(distVector), "all") - 1;
                for m = 1:numNodes % 遍历所有相邻路由器（接收到路由信息，触发路由更新）
                    if isinf(G(m, n)) % 不相邻
                        continue
                    end
                    transferCount = transferCount + 1;
                    transferSize = transferSize + dvSize;
                end
            end
        end

        % 路由收敛
        if ~updated
            fprintf('路由收敛，总传输次数 %d 次，总传输包大小 %d 。\n', transferCount, transferSize);
            break
        end

        loopCount = loopCount + 1;
        % 输出路由表
        fprintf('---------------- LOOP %d ----------------\n', numLoop)
        disp('距离矩阵:');
        disp(distMatrix);
        disp('下一跳矩阵:')
        disp(nextHopMatrix);
    end
    fprintf('达到循环次数上限或路由收敛，路由信息计算结束。\n');
end