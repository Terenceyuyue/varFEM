function output = int1dvec(Th,Coef,Trial,Test,feSpace,quadOrder)

if nargin == 4, feSpace = 'P1'; quadOrder = 3; end % default: P1
if nargin == 5, quadOrder = 3; end

%% -------------------------- Bilinear form --------------------------

%% -------------------------- Linear form --------------------------
if isa(Coef,'function_handle')
    trf = eye(2); f = Coef;
    Coef1 = @(pz) f(pz)*trf(:, 1);   Coef2 = @(pz) f(pz)*trf(:, 2); 
    F1 = int1d(Th,Coef1,[],Test,feSpace,quadOrder);     
    F2 = int1d(Th,Coef2,[],Test,feSpace,quadOrder);
else
    Coef1 = Coef{1}; Coef2 = Coef{2}; % matrix
    F1 = int1d(Th,Coef1,[],Test,feSpace,quadOrder);
    F2 = int1d(Th,Coef2,[],Test,feSpace,quadOrder);
end
output = [F1; F2];



