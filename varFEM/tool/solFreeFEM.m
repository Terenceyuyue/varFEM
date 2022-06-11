function u = solFreeFEM(filename)

u = importdata(filename);
u = u'; u = u(:);
ind = ~isnan(u);
u = u(ind);
u(1) = [];