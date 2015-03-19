
% $Id$

function [] = GaussLegendreTest()

n = 199; 

[ x, w ] = GaussLegendre(n);

f = myf(x);

plot(x, real(f), x, imag(f), 'LineWidth', 4);
hold on;

m = 100;

P = LegendreP2(m, x);

L = [0:1:m] + 0.5;

fl = (w.*f)'*P.*L;

f2 = P*fl';

max(abs(real(f-f2)))
max(abs(imag(f-f2)))

plot(x,real(f2), 'ro-', x,imag(f2), 'ro-');
hold off;

return

function [ v ] = myf(x) 
v = 1./sqrt(3+x-x.*x).*exp(-cos(6*x)+2*x.^2 - x).*sin(10*x).*exp(25*j*x);
return
