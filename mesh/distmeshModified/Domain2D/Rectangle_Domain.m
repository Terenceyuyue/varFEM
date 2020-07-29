function [fd,fh,BdBox,pfix] = Rectangle_Domain

a = 0; b = 1; c = 0; d = 2;

fd = @(p) drectangle(p,a,b,c,d);
fh = @huniform;
BdBox = [a b c d];
pfix = [a c; b c; b d; a d];
