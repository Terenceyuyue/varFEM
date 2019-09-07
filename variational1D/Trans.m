function cw = Trans(w,xa,xb)
    % coefficients of the affine transformation
    if strcmpi(w,'u.val') || strcmpi(w,'v.val')  
        nel = length(xa); cw = ones(nel,1);
    elseif strcmpi(w,'u.dx') || strcmpi(w,'v.dx')
        cw = 1./(xb-xa);
    else
        error("The input character must be 'u.val', 'u.dx', 'v.val' or 'v.dx'.");
    end
