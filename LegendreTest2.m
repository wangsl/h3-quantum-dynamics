
% $Id$

function LegendreTest2
%clear all
%clc
%close all

%format long

global x
global w
global F1
global F2
global G
global a
global P
global n

% for discrete data, we have to use grids from Legendre polynomial roots

n = 190;

[ x, w ] = GaussLegendre(n);
F1 = myf(x);

subplot(2,1,1);
plot(x, real(F1), x, imag(F1), 'LineWidth', 3);
hold on;

tic
P = LegendreP2(n, x);
a = alpha(x, w, F1);
toc

tic
F2 = zeros(size(x));
for L = 1 : n
  F2 = F2 + a(L)*P(:,L);
end
toc

subplot(2,1,1);
plot(x, real(F2), 'r--', x, imag(F2), 'r--', 'LineWidth', 2);
hold off;

tic
G = zeros(size(x));
for L = 1 : n
  G = G - (L-1)*L*a(L)*P(:,L);
end
toc

subplot(2,1,2);
plot(x, real(G), x, imag(G), 'LineWidth', 2);
hold off

print -dpdf GL-exp-test-2.pdf

min(abs(F1-F2))

return

%%%

function [ v ] = myf(x) 
v = 1./sqrt(3+x-x.*x).*exp(-cos(6*x)+2*x.^2 - x).*sin(10*x).* ...
    exp(25*j*x);
return

%%%

function [ v ] = alpha(x, w, F1)

global P;

n = numel(x);
assert(n == numel(w));
assert(n == numel(F1));

v = zeros(n, 1);
for L = 1 : n 
  v(L) = sum(w.*F1.*P(:,L));
  v(L) = (L-1+1/2)*v(L);
end

return



