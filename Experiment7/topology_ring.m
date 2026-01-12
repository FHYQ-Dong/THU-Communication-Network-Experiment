clc
clear
% 环状拓扑构建
numNodes = 8;
weightMax = 10;
s = 1:numNodes;
t = circshift(s,-1);
w = randi(weightMax,1,numNodes);
%  节点编号   1    2    3    4    5    6    7    8
nodeNames = {'a', 'b', 'c', 'd', 'e', 'f', 'g','h'};
G = graph(s, t, w, nodeNames);
% 固定可视化时节点位置
nodePositions = zeros(numNodes,2);
deltaTheta = 2*pi/numNodes;
for k = 1:numNodes
nodePositions(k,1) = real(exp(1i*(k-1)*deltaTheta));
nodePositions(k,2) = imag(exp(1i*(k-1)*deltaTheta));
end
% 图拓扑可视化
figure;
plot(G, 'XData', nodePositions(:, 1), 'YData', nodePositions(:, 2), 'EdgeLabel', G.Edges.Weight);
title('节点拓扑');
axis equal;
hold on;
save('topology_ring.mat', 's', 't', 'w', 'nodePositions', 'nodeNames');