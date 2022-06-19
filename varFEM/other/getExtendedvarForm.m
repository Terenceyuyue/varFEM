function [Coef,Test,Trial] = getExtendedvarForm(Coef,Test,Trial)
%% GETEXTENDEDVARFORM expands the blinear forms
%
% e.g. Eij(u)Eij(v) = v1x*u1x + v2x*u2y + 0.5(v1y+v2x)*(u1y+u2x)
% is expanded as
%  Coef  = {1, 1, 0.5, 0.5, 0.5, 0.5};
%  Test  = {'v1.dx', 'v2.dy', 'v1.dy', 'v1.dy', 'v2.dx', 'v2.dx'};
%  Trial = {'u1.dx', 'u2.dy', 'u1.dy', 'u2.dx', 'u1.dy', 'u2.dx'};
%

% check the existence of linear combination
idv = cellfun(@(str) mycontains(str,'+'), Test);  idv = sum(idv);
idu = cellfun(@(str) mycontains(str,'+'), Trial); idu = sum(idu);
if ~idv && ~idu, return; end  % no '+' in Test and Trial

cc = cell(1,length(Coef)); vv = cc; uu = cc;
for ss = 1:length(Test)
    strv = Test{ss}; stru = Trial{ss};
    [cs, vs, us] = getvarForm(strv, stru);
    cf = Coef{ss};
    for i = 1:length(cs)
        if isnumeric(cf), cs{i} = cf*cs{i}; end
        if isa(cf, 'function_handle'), cs{i} = @(p) cf(p)*cs{i}; end
    end
    cc{ss} = cs; vv{ss} = vs; uu{ss} = us;
end

Coef  = horzcat(cc{:});
Test  = horzcat(vv{:});
Trial = horzcat(uu{:});
end


function [Coef,Test,Trial] = getvarForm(strv, stru)
% example: strv = 'v1.dx + 2*v2.val';
%          stru = '3*u1.val + pde.a*u2.dx + u2.val';

% remove spaces
strv = strrep(strv,' ','');  stru = strrep(stru,' ','');
% pairs: split on '+'
strv = strsplit(strv, '+');
stru = strsplit(stru, '+');
% triples: split on '*'
nv    = length(strv);    nu    = length(stru);
cvv   = cell(1,nv);      cuu   = cell(1,nu);
strvv = cell(1,nv);      struu = cell(1,nu);
for i = 1:nv
    str = strv{i}; str = strsplit(str,'*');
    if length(str)==1, cvv{i} = 1;                  strvv{i} = str{1}; end
    if length(str)==2, cvv{i} = str2double(str{1}); strvv{i} = str{2}; end
end
for i = 1:nu
    str = stru{i}; str = strsplit(str,'*');
    if length(str)==1, cuu{i} = 1;                  struu{i} = str{1}; end
    if length(str)==2, cuu{i} = str2double(str{1}); struu{i} = str{2}; end
end
% cross pairs: Coef, Test, Trial
[Coef,Test,Trial] = deal(cell(1,nu*nv)); 
s = 1;
for i = 1:nu
    for j = 1:nv
        Coef{s} = cvv{i}*cuu{j};
        Test{s} = strvv{i};
        Trial{s} = struu{j};
        s = s + 1;
    end
end
end
