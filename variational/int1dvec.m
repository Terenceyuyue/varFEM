function output = int1dvec(Th,Coef,Trial,Test,feSpace,quadOrder)

if nargin == 4, feSpace = 'P1'; quadOrder = 3; end % default: P1
if nargin == 5, quadOrder = 3; end

%% -------------------------- Bilinear form --------------------------
if ~isempty(Trial)
    return; % the remaining code will be neglected.
end

%% -------------------------- Linear form --------------------------

if strcmpi(Test, 'v1.val') % v1.val
    F1 = int1d(Th,Coef,[],Test,feSpace,quadOrder);
    NNdof = size(F1,1);  F = zeros(2*NNdof,1);
    F(1:NNdof) = F(1:NNdof) + F1;
end
if strcmpi(Test, 'v2.val') % v2.val
    F2 = int1d(Th,Coef,[],Test,feSpace,quadOrder);
    NNdof = size(F2,1);  F = zeros(2*NNdof,1);
    F(NNdof+1:end) = F(NNdof+1:end) + F2;
end
if strcmpi(Test, 'v.val')  % v.val = [v1.val, v2.val]
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
    F = [F1; F2];
end
output = F;