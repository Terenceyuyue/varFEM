function [w1,w2] = Refbase1D(w)
    % Trial or Test functions on reference interval [0,1]
    if strcmpi(w,'u.val') || strcmpi(w,'v.val')
        w1 = @(x) 1-x;  w2 = @(x) x;
    elseif strcmpi(w,'u.dx') || strcmpi(w,'v.dx')
        w1 = @(x) -ones(size(x));  w2 = @(x) ones(size(x));
    else
        error("The input character must be 'u.val', 'u.dx', 'v.val' or 'v.dx'.");
    end