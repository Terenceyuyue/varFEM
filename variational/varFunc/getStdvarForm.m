function [Test,Trial] = getStdvarForm(vstr, ustr, Test, Trial)
%% GETSTDVARFORM transforms the user-defined symbols to the standard ones
%
%  vstr --> [v1,v2,q] --> [v1,v2,v3]
%  stru --> [u1,u2,p] --> [u1,u2,u3]
% e.g.
% vstr = {'v1','v2','q'}; ustr = {'u1','u2','p'};
% Test  = {'v1.grad', 'v2.grad', 'v1.dx',  'v2.dy',  'q.val', 'q.val', 'q.val'};
% Trial = {'u1.grad', 'u2.grad', 'p.val',  'p.val', 'u1.dx',  'u2.dy',  'p.val'};
%

for i = 1:length(vstr)
    vi = sprintf('v%d',i);  ui = sprintf('u%d',i);
    for s = 1:length(Test)
        strv = Test{s}; stru = Trial{s};
        Test{s} = strrep(strv, vstr{i},  vi);
        Trial{s} = strrep(stru, ustr{i}, ui);
    end
end

