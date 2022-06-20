function [node,elem] = getMeshFreeFEM(meshName)

fid = fopen(meshName,'r');

%% Numbers 
Num = fscanf(fid, '%d', 3);  % [N, NT, NbdEdge]
N = Num(1); NT = Num(2); 

%% node
node = fscanf(fid, '%f%f', [3,N]);
node = node(1:2,:)';

%% elem
elem = fscanf(fid, '%d%d%d', [4,NT]);
elem = elem(1:3,:)';

fclose(fid);