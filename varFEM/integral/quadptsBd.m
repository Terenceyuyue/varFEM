function [lambdaBd,weightBd] = quadptsBd(order)
%%Gauss quadrature points along the boundary of triangles

% quadrature on [0,1]
[lambda1d,weight1d] = quadpts1(order);  ng = length(weight1d);
[~,id] = sort(lambda1d(:,1)); 
lambda1d = lambda1d(id,:); weight1d = weight1d(id);

lambdae1 = [zeros(ng,1), lambda1d(:,2), lambda1d(:,1)];
lambdae2 = [lambda1d(:,1), zeros(ng,1), lambda1d(:,2)];
lambdae3 = 1 - lambdae1 - lambdae2;
lambdaBd = [lambdae1; lambdae2; lambdae3];
weightBd = repmat(weight1d,1,3);