clc
clear
% 网状拓扑构建
numNodes = 8;
weightMax = 10;
numEdges = numNodes*(numNodes-1)/2;
s = zeros(1,numEdges);
lastIdx = 1;
for i = 1:numNodes-1
s(lastIdx:lastIdx+numNodes-i-1) = i;
lastIdx = lastIdx + numNodes - i;
end
t = zeros(1,numEdges);
lastIdx = 1;
for i = 1:numNodes-1
for j = i+1:numNodes
t(lastIdx) = j;
lastIdx = lastIdx + 1;
end
end
rng(1);
w = randi(weightMax,1,numEdges);
%  节点编号   1    2    3    4    5    6    7    8
nodeNames = {'a', 'b', 'c', 'd', 'e', 'f', 'g','h'};
G = graph(s, t, w, nodeNames);
% 固定可视化时节点位置
nodePositions = [0 0; 0 1; 0.5 1.5; 2 0; 2 2; 4 0.5; 4 1; 4 1.75];
% 图拓扑可视化
figure;
plot(G, 'XData', nodePositions(:, 1), 'YData', nodePositions(:, 2), 'EdgeLabel', G.Edges.Weight);
xlim([-0.2;4.2]);
ylim([-0.2,2.2]);
title('节点拓扑');
hold on;
save('topology_mesh.mat', 's', 't', 'w', 'nodePositions', 'nodeNames');