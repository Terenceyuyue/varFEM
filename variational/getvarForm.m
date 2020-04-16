function [Coef,Trial,Test] = getvarForm(stru, strv)
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
    % cross pairs: Coef, Trial, Test
    % Coef
    cvv = repmat(cvv, 1, nu);
    Coef1 = cvv([1:2:end, 2:2:end]);
    Coef2 = repmat(cuu, 1, nv);
    Coef = cellfun(@(x,y) x*y, Coef1, Coef2, 'UniformOutput', false);
    % Test
    strvv = repmat(strvv, 1, nu);
    Test = strvv([1:2:end, 2:2:end]);
    % Trial
    Trial = repmat(struu, 1, nv);
