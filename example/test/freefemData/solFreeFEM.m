function u = solFreeFEM(filename)

fid = fopen(filename,'r');
N = fscanf(fid, '%d', 1);
u = fscanf(fid, '%f', N);
fclose(fid);