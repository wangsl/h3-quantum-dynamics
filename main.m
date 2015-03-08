
% $Id$
%
clear all
%close all

%clc
format long

global UseLSTH
%global psi

setenv('OMP_NUM_THREADS', '20');

% UseLSTH = true;

MassAU = 1.822888484929367e+03;

mH = 1.007825*MassAU;
mD = 2.01410178*MassAU;
mT = 3.0160492*MassAU;

masses = [ mH mH mH ];

vH2Min = -0.174495770896975;

r1.n = int32(512);
r1.r = linspace(0.4, 14.0, r1.n);
r1.dr = r1.r(2) - r1.r(1);
r1.mass = 2*mH/3;
r1.r0 = 10.0;
r1.k0 = 4.0;
r1.delta = 0.12;

r2.n = int32(256);
r2.r = linspace(0.4, 10.0, r2.n);
r2.dr = r2.r(2)-r2.r(1);
r2.mass = mH/2;

theta.n = int32(180);
theta.m = int32(120);
[ theta.x, theta.w ] = GaussLegendre(theta.n);

pot = H3PESJacobi(r1.r, r2.r, acos(theta.x), masses);

jRot = 6;
nVib = 4;
%psi = InitWavePacket(r1, r2, theta, jRot, nVib);

%psi = zeros(4,2,2)
%psi(1,1,1) = 1;
%psi(2,1,1) = 2;

%psi

psi = InitWavePacket(r1, r2, theta, jRot, nVib);

tic
TimeEvolutionMex(r1, r2, theta, pot, psi)
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


