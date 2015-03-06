
% $Id$

clear all
%close all

clc
format long

global x
global F1
global F2
global G

LegendreTest

for i = 1 : 5 : numel(x)
  fprintf(1, '%12.8f %12.8f %12.6e %12.8f %12.8f\n', ...
	  x(i), F1(i), abs(F1(i)-F2(i)), real(G(i)), imag(G(i)));
end


