

function BooleTest

clear all
clc
format long

J = 10;

n = 20*(J-1)+1
x = linspace(-1, 1, n);
dx = x(2) - x(1)

L = 20;
P = LegendreP(L, x);

plot(x, P);

P = P.^2;

tic
s = Boole(P(1:J), dx, J);;
i = J;
while true 
  j = i+J-1;
  s = s + Boole(P(i:j), dx, J);
  i = j;
  if j >= n
    break
  end
end
s = s*(L+1/2)
toc

tic
sum(P)*dx*(L+1/2)
toc
