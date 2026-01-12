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

% 迭代更新距离矩阵和下一跳矩阵
[distMatrix, nextHopMatrix] = updateRoute(G, distMatrix, nextHopMatrix, 100);

% 链路失效后路由收敛过程与"Count to infinity"问题
bads = 6;% 失效链路源路由节点
badt = 7;% 失效链路另一路由节点
fprintf('\n%d-%d链路断开...\n\n', bads, badt);
G(bads, badt) = inf;
G(badt, bads) = inf;
% 清空两个路由节点路由表中通过对方的路由信息（距离矢量和下一跳）
passBadT = (nextHopMatrix(bads, :) == badt);
passBadS = (nextHopMatrix(badt, :) == bads);
distMatrix(bads, passBadT) = inf;
distMatrix(badt, passBadS) = inf;
nextHopMatrix(bads, passBadT) = 0;
nextHopMatrix(badt, passBadS) = 0;

fprintf('重新计算路由...\n');
[distMatrix, nextHopMatrix] = updateRoute(G, distMatrix, nextHopMatrix, 16);


% 距离向量法（Bellman-Ford 算法）
% 根据拓扑G，更新路由表（距离矩阵和下一跳矩阵）直到收敛
% 返回收敛的路由表（距离矩阵和下一跳矩阵）
function [distMatrix, nextHopMatrix] = updateRoute(G, distMatrix, nextHopMatrix, loops)
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

    % 迭代更新距离矩阵和下一跳矩阵
    for numLoop = 1:loops % 不断迭代更新直到收敛
        updated = false; % 是否收敛
        updateID = zeros(1, numNodes); % 节点在一次循环后是否更新路由表
        % 假设每个路由器定时向相邻路由器广播路由信息，广播的顺序与路由器节点编号顺序一致
        for n = 1:numNodes % 遍历所有节点（进行广播）
            % 从n发出的距离矢量
            distVector = distMatrix(n, :);
            % 本次循环以下代码，模拟节点n的本轮广播中向所有其它节点m传输的开销统计以及其它节点m在收到节点n广播本轮距离矢量之后的动作。
            % 仿真的其它节点动作，包含这些节点在内部对距离矢量和路由表的更新，不包含对外广播距离矢量。
             
            % 距离矢量表大小（inf的视为不存在）
            dvSize = sum(~isinf(distVector), "all") - 1;
            for m = 1:numNodes % 遍历所有相邻路由器（接收到路由信息，触发路由更新）
                if isinf(G(m, n)) % 不相邻
                    continue
                end

                % 第一次统计：传输次数+=1，传输包大小+=距离矢量个数
                if numLoop == 1
                    transferCount = transferCount + 1;
                    transferSize = transferSize + dvSize;
                end

                % 路由表更新：距离矢量从n传输到m，m接收后进行出力
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
            % 若节点的路由表更新，应该立刻广播给邻居节点
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
