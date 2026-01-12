clc
clear
% 星状拓扑构建
numNodes = 9;
weightMax = 10;
s = ones(1,numNodes-1);
t = 2:numNodes;
rng(1);
w = randi(weightMax,1,numNodes-1);
%  节点编号   1    2    3    4    5    6    7    8   9
nodeNames = {'a', 'b', 'c', 'd', 'e', 'f', 'g','h','i'};
G = graph(s, t, w, nodeNames);
% 固定可视化时节点位置
nodePositions = [1 1; 0 1; 0 2; 1 2; 2 2; 2 1; 2 0; 1 0; 0 0];
% 图拓扑可视化
figure;
plot(G, 'XData', nodePositions(:, 1), 'YData', nodePositions(:, 2), 'EdgeLabel', G.Edges.Weight);
title('节点拓扑');
axis equal;
hold on;
save('topology_star.mat', 's', 't', 'w', 'nodePositions', 'nodeNames');