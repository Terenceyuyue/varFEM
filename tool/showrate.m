function r = showrate(h,err,opt1,opt2)

err(err == 0) = 1e-16; % Prevent the case err = 0, log(err) = -Inf.
p = polyfit(log(h(1:end)),log(err(1:end)),1);
r = p(1);
s = 0.75*err(1)/h(1)^r;

loglog(h,err,opt1,'linewidth',2);
hold on
loglog(h,s*h.^r,opt2,'linewidth',1);