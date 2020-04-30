function showrateh(h,ErrL2,ErrH1,str1,str2)
if nargin==3
    str1 = '||u - u_h||';   % L2
    str2 = '||Du - Du_h||'; % H1    
end

r2 = showrate(h,ErrH1,'r-*','k.-');
hold on
r1 = showrate(h,ErrL2,'b-s','k--');

h_legend = legend(str2,['O (h^{' num2str(r2,2) '})'],...
    str1,['O (h^{' num2str(r1,2) '})'],'location','best');
set(h_legend,'FontSize',10);