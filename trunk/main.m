
% $Id$

clear all
%close all

clc
format long

global masses
global UseLSTH

setenv('OMP_NUM_THREADS', '20');

% UseLSTH = true;

%ElectronMass = 9.10938291E-31;
%AMU = 1.66053892E-27;
MassAU = 1.822888484929367e+03;

mH = 1.007825*MassAU;
mD = 2.01410178*MassAU;
mT = 3.0160492*MassAU;

masses = [ mH mH mH ];

vH2Min = -0.174495770896975;

r.n = 512;
r.r = linspace(0.4, 16.0, r.n);
r.dr = r.r(2) - r.r(1);
r.mass = mH/2;

energy_levels = [ 0 1 10 13 14 15 16 ];

[ e, psi ] = H2WaveFunction(r, energy_levels);

plot(r.r, psi, 'LineWidth', 2);
legend(arrayfun(@num2str, energy_levels, 'UniformOutput', false));

print -dpdf H2-wavefunction.pdf

return



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


