
% $Id$

function [] = GaussLegendreTest()

n = 190;

[ x, w ] = GaussLegendre(n);

f = myf(x);

plot(x, real(f), x, imag(f), 'ro-');

for L = 1 : 60
  P = LegendreP(L-1,x);
  s = sum(w.*f.*P);
  fprintf(1, '%4d %20.15f%20.15f\n', L-1, real(s), imag(s));
end

return

function [ v ] = myf(x) 
v = 1./sqrt(3+x-x.*x).*exp(-cos(6*x)+2*x.^2 - x).*sin(10*x).*exp(25*j*x);
return
