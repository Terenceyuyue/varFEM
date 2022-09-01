function showrateErr(h,varargin)
%showrateErr displays convergence rates of an err sequence with given names
% such as ErrL2, ErrH1, ErrH2, ErrI
%  e.g. showrateErr(h,ErrL2,ErrH1,ErrH2,ErrI)
%

n = nargin-1;

%% Define the line styles
lab1 = {'r-*', 'b-s', 'b-o', 'k-d'};
lab2 = {'k.-', 'k--', 'r.-', 'k:.'};

%% Plot convergence rates according to the inputnames
str = cell(1,2*n);
for i = 1:n
    stri = inputname(i+1);
    r = showrate(h,varargin{i},lab1{i},lab2{i});
    str{2*i-1} = stri;
    str{2*i} = ['O (h^{',  num2str(r,2),  '})'];
end

h_legend = legend(str,'location','best');
set(h_legend,'FontSize',10);