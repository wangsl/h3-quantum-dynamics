
% $Id$

clear all
%close all

clc
format long

global masses

ElectronMass = 9.10938291E-31;
AMU = 1.66053892E-27;
MassAU = AMU/ElectronMass;

mH = 1.007825*MassAU;
mD = 2.01410178*MassAU;
mT = 3.0160492*MassAU;

masses = [ mH mH mH ];

R = linspace(0.2, 14.0, 256);
r = linspace(0.2, 14.0, 256);
Theta = 90.0/180.0*pi;

tic
V = H3JacobiPES(R, r, Theta);
toc

[ ~, hPES ] = contourf(R, r, V', [ -0.1:0.01:-0.001 0.001:0.01:0.3 ]);
set(hPES, 'LineWidth', 0.75);
set(hPES, 'LineColor', 'black');
colorbar('vert');
axis square

return
