
% $Id$

clear all
%close all
clc

format long

global UseLSTH

setenv('OMP_NUM_THREADS', '20');

% UseLSTH = true;

MassAU = 1.822888484929367e+03;

mH = 1.007825*MassAU;
mD = 2.01410178*MassAU;
mT = 3.0160492*MassAU;

masses = [ mH mH mH ];

vH2Min = -0.174495770896975;

% time

time.total_steps = int32(800);
time.time_step = 1.0;
time.steps = int32(0);

% r1: R

r1.n = int32(128);
r1.r = linspace(0.3, 14.0, r1.n);
r1.dr = r1.r(2) - r1.r(1);
r1.mass = 2*mH/3;
r1.r0 = 4.0;
r1.k0 = 1.0;
r1.delta = 0.12;

dump1.Cd = 3.0;
dump1.xd = 12.0;
dump1.dump = WoodsSaxon(dump1.Cd, dump1.xd, r1.r);

% r2: r

r2.n = int32(128);
r2.r = linspace(0.3, 12.0, r2.n);
r2.dr = r2.r(2) - r2.r(1);
r2.mass = mH/2;

dump2.Cd = 3.0;
dump2.xd = 10.0;
dump2.dump = WoodsSaxon(dump2.Cd, dump2.xd, r2.r);

% angle:

theta.n = int32(80);
theta.m = int32(60);
[ theta.x, theta.w ] = GaussLegendre(theta.n);

theta.legendre = LegendreP2(double(theta.m), theta.x);
theta.legendre = theta.legendre';

% options

options.wave_to_matlab = 'C2Matlab.m'

pot = H3PESJacobi(r1.r, r2.r, acos(theta.x), masses);

jRot = 2;
nVib = 4;

[ psi, eH2, psiH2 ] = InitWavePacket(r1, r2, theta, jRot, nVib);

%H2WaveFunctionAnalysis(r2.r, psiH2, r2.mass, jRot)
%eKGaussian = 1/(2*r1.mass)*(r1.k0^2 + 1/(2*r1.delta^2))

%PlotPotWave(r1, r2, pot, psi)

tic
TimeEvolutionMex(r1, r2, theta, pot, psi, time, options, dump1, dump2)
toc

tic
PSI = psi(1:2:end, :, :) + j*psi(2:2:end, :, :);
a = sum(sum(conj(PSI).*PSI));
a = reshape(a, [numel(a), 1]);
sum(theta.w.*a)*r1.dr*r2.dr
toc

return

pot = H3PESJacobi(r1.r, r2.r, acos(theta.x), masses);

jRot = 6;
nVib = 4;
psi = InitWavePacket(r1, r2, theta, jRot, nVib);

PSI = psi(1:2:end, :, :) + j*psi(2:2:end, :, :);

a = sum(sum(conj(PSI).*PSI));
a = reshape(a, [numel(a), 1]);

sum(theta.w.*a)*r1.dr*r2.dr

clearvars PSI a

return

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


