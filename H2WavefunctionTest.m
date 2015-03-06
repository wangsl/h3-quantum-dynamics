
% $Id$

clear all
%close all

clc
format long

global masses
global UseLSTH

setenv('OMP_NUM_THREADS', '20');

% UseLSTH = true;

MassAU = 1.822888484929367e+03;

mH = 1.007825*MassAU;
mD = 2.01410178*MassAU;
mT = 3.0160492*MassAU;

masses = [ mH mH mH ];

vH2Min = -0.174495770896975;

r.n = 4096;
r.r = linspace(0.4, 16.0, r.n);
r.dr = r.r(2) - r.r(1);
r.mass = mH/2;

energy_levels = [ 0 1 6 13 14 15 16 ];

[ e, psi ] = H2Wavefunction(r, energy_levels);

plot(r.r, psi, 'LineWidth', 2);
legend(arrayfun(@num2str, energy_levels, 'UniformOutput', false));

print -dpdf H2-wavefunction.pdf

return
