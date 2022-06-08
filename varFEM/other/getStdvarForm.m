function varargout = getStdvarForm(varargin)
%function [Test,Trial] = getStdvarForm(vstr, Test, ustr, Trial)
%function Test = getStdvarForm(vstr, Test)
%%getStdvarForm transforms the user-defined symbols to the standard ones
%
%  vstr --> [v1,v2,q] --> [v1,v2,v3]
%  stru --> [u1,u2,p] --> [u1,u2,u3]
% e.g.
% vstr = {'v1','v2','q'}; ustr = {'u1','u2','p'};
% Test  = {'v1.grad', 'v2.grad', 'v1.dx',  'v2.dy',  'q.val', 'q.val', 'q.val'};
% Trial = {'u1.grad', 'u2.grad', 'p.val',  'p.val', 'u1.dx',  'u2.dy',  'p.val'};
%

nout = nargin/2;
for s = 1:nout
    str = varargin{2*s-1};  seq = varargin{2*s};
    if mycontains(lower(inputname(2*s)), 'test')
        r = 'v';
    else
        r = 'u';
    end
    for i = 1:length(str)
        vi = [r,num2str(i)];  
        for ss = 1:length(seq)
            strv = seq{ss}; 
            seq{ss} = strrep(strv, str{i},  vi);
        end
    end
    varargout{s} = seq; %#ok<AGROW> 
end