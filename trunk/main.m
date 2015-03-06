
% $Id$

clear all
%close all

clc
format long

global masses
global UseLSTH

setenv('OMP_NUM_THREADS', '20');

% UseLSTH = true;

ElectronMass = 9.10938291E-31;
AMU = 1.66053892E-27;
MassAU = AMU/ElectronMass;

mH = 1.007825*MassAU;
mD = 2.01410178*MassAU;
mT = 3.0160492*MassAU;

masses = [ mH mH mH ];

r.n = 4096;
r.r = linspace(0.4, 20.0, r.n);
r.dr = r.r(2) - r.r(1);
r.mass = mH/2;

[ e, psi ] = H2Wavefunction(r, 0);

return


R = linspace(0.45, 14.0, 512);
r = linspace(0.45, 8.0, 512);
Theta = linspace(0.0, pi, 180);

tic
V = H3PESJacobi(R, r, Theta);
toc

V1 = V(:, :, 1);

[ ~, hPES ] = contour(R, r, V1', [ -0.1:0.01:-0.001 0.001:0.01:0.4 ]);
set(hPES, 'LineWidth', 0.75);
set(hPES, 'LineColor', 'black');
colorbar('vert');
%axis square
%hold on;


