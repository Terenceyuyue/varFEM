function [fd,fh,BdBox,pfix] = Lshape_Domain

a = 0; b = 1; c = 0; d = 1;
xm = (a+b)/2;  ym = (c+d)/2;

fd = @fd1;
fh = @huniform;
BdBox = [a b c d];
pfix = [a c; xm c; xm ym; b ym; b d; a d];

    function dist = fd1(p)        
        d1 = drectangle(p,a,b,c,d);
        d2 = drectangle(p,xm,b,c,ym);
        dist = ddiff(d1,d2);
    end
end