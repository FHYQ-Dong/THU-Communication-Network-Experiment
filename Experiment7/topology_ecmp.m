clc
clear
% 菱形拓扑构建
s = [1 1 2 3];
t = [2 3 4 4];
w = [1 1 1 1];
%  节点编号   1    2    3    4
nodeNames = {'a', 'b', 'c', 'd'};
G = graph(s, t, w, nodeNames);
% 固定可视化时节点位置
nodePositions = [0 1; 1 2; 1 0; 2 1];
% 图拓扑可视化
figure;
plot(G, 'XData', nodePositions(:, 1), 'YData', nodePositions(:, 2), 'EdgeLabel', G.Edges.Weight);
title('等价多路径(ECMP)菱形拓扑');
axis equal;
hold on;
save('topology_ecmp.mat', 's', 't', 'w', 'nodePositions', 'nodeNames');