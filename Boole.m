
% $Id$

% Ref: Abramowitz and Stegun (1964) P886 Eq. 25.4.14 to 25.4.20

function [ s ] = Boole(f, dx, n)
switch n 
 case 5
  s = Boole5(f, dx);
 case 6
  s = Boole6(f, dx);
 case 7
  s = Boole7(f, dx);
 case 8
  s = Boole8(f, dx);
 case 9
  s = Boole9(f, dx);
 case 10
  s = Boole10(f, dx);
 otherwise
  error('no support')
end

return

% 

function [ s ] = Boole5(f, dx)
s = 7*(f(1)+f(5)) + 32*(f(2)+f(4)) + 12*f(3);
s = s*2*dx/45;

function [ s ] = Boole6(f, dx)
s = 19*(f(1)+f(6)) + 75*(f(2)+f(5)) + 50*(f(3)+f(4));
s = s*5*dx/288;

function [ s ] = Boole7(f, dx)
s = 41*(f(1)+f(7)) + 216*(f(2)+f(6)) + 27*(f(3)+f(5)) + 272*f(4);
s = s*dx/140;

function [ s ] = Boole8(f, dx)
s = 751*(f(1)+f(8)) + 3577*(f(2)+f(7)) + 1323*(f(3)+f(6)) + ...
    2989*(f(4)+f(5));
s = s*7*dx/17280;

function [ s ] = Boole9(f, dx)
s = 989*(f(1)+f(9)) + 5888*(f(2)+f(8)) - 928*(f(3)+f(7)) + ...
    10496*(f(4)+f(6)) - 4540*f(5);
s = s*4*dx/14175;

function [ s ] = Boole10(f, dx)
s = 2857*(f(1)+f(10)) + 15741*(f(2)+f(9)) + 1080*(f(3)+f(8)) + ...
    19344*(f(4)+f(7)) + 5778*(f(5)+f(6));
s = s*9*dx/89600;
