function g = ellipseg(a,b)

g = [4   4   4   4  % 4 for ellipse
    a   0  -a   0   % x: starting points
    0  -a   0   a   % x: ending points
    0   b   0  -b   % y: starting points
    b   0  -b   0   % y: ending points
    1   1   1   1   % label of subdomain on the left
    0   0   0   0   % label of subdomain on the right
    0   0   0   0   % center x
    0   0   0   0   % center y
    a   a   a   a   % semimajor axis a
    b   b   b   b   % semiminor axis b
    0   0   0   0   % rotation angle
    ];