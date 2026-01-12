clc
clear
% 图拓扑构建
s = [1 1 1 2 2 3 3 3 4 5 6];
t = [2 3 4 3 4 4 5 6 5 6 7];
w = [2 5 1 3 2 3 1 5 1 2 1];
%  节点编号   1    2    3    4    5    6    7 
nodeNames = {'a', 'b', 'c', 'd', 'e', 'f', 'g'};
G = graph(s, t, w, nodeNames);
% 固定可视化时节点位置
nodePositions = [0 1; 1 1; 2 2; 1 0; 2 0; 3 1; 4 1];
% 图拓扑可视化
figure;
plot(G, 'XData', nodePositions(:, 1), 'YData', nodePositions(:, 2), 'EdgeLabel', G.Edges.Weight);
title('节点拓扑');
axis equal;
hold on;
save('topology_general.mat', 's', 't', 'w', 'nodePositions', 'nodeNames');