  
% $Id$

clear all
%close all
clc

format long

global UseLSTH
global H2eV 

setenv('OMP_NUM_THREADS', '20');

% UseLSTH = true;

%MassAU = 1.822888484929367e+03;

ElectronMass = 9.10938291E-31;
AMU = 1.66053892E-27;
MassAU = AMU/ElectronMass;

mH = 1.007825*MassAU;
mD = 2.01410178*MassAU;
mT = 3.0160492*MassAU;

masses = [ mH mH mH ];

vH2Min = -0.174495770896975;

H2eV = 27.21138505;

% time

time.total_steps = int32(4000);
time.time_step = 5;
time.steps = int32(0);

% r1: R

r1.n = int32(256);
r1.r = linspace(0.4, 14.0, r1.n);
r1.dr = r1.r(2) - r1.r(1);
r1.mass = 2*mH/3;
r1.r0 = 8.0;
r1.k0 = 5.0;
r1.delta = 0.2;

dump1.Cd = 3.0;
dump1.xd = 12.0;
dump1.dump = WoodsSaxon(dump1.Cd, dump1.xd, r1.r);

% r2: r

r2.n = int32(256);
r2.r = linspace(0.4, 12.0, r2.n);
r2.dr = r2.r(2) - r2.r(1);
r2.mass = mH/2;

% dump functions

dump2.Cd = 3.0;
dump2.xd = 10.0;
dump2.dump = WoodsSaxon(dump2.Cd, dump2.xd, r2.r);

% dividing surface

divSurf.rd = 7.0;
divSurf.n = int32((divSurf.rd - min(r2.r))/r2.dr);
r2Div = double(divSurf.n)*r2.dr + min(r2.r);
fprintf(' Dviding surface: %.8f\n', r2Div);

% angle:

dimensions = 3;

if dimensions == 2 
  % for 2 dimensional case
  theta.n = int32(1);
  theta.m = int32(0);
  theta.x = 1.0;
  theta.w = 2.0;
else 
  % for 3 dimensional case
  theta.n = int32(120);
  theta.m = int32(110);
  [ theta.x, theta.w ] = GaussLegendre(theta.n);
end
  
theta.legendre = LegendreP2(double(theta.m), theta.x);
% transpose Legendre polynomials in order to do 
% matrix multiplication in C++ and Fortran LegTransform.F
theta.legendre = theta.legendre';

% options

options.wave_to_matlab = 'C2Matlab.m';

pot = H3PESJacobi(r1.r, r2.r, acos(theta.x), masses);

jRot = 0;
nVib = 1;

[ psi, eH2, psiH2 ] = InitWavePacket(r1, r2, theta, jRot, nVib);

% cummulative reaction probabilities

CRP.eH2 = eH2;
CRP.n_gradient_points = 11;
CRP.eLeft = 0.4/H2eV - eH2;
CRP.eRight = 4.0/H2eV - eH2;
CRP.energy = linspace(CRP.eLeft, CRP.eRight, 400);
CRP.etaSq = EtaSq(r1, CRP.energy);
CRP.faiE = complex(zeros(r1.n, theta.n, length(CRP.energy)));
CRP.DfaiE = complex(zeros(r1.n, theta.n, length(CRP.energy)));
CRP.CRP = zeros(size(CRP.energy));

% wrapper the data will not used in C++ as others
 
others.divSurf = divSurf;
others.CRP = CRP;

% time evolution

tic
TimeEvolutionMex(r1, r2, theta, pot, psi, time, options, ...
		 dump1, dump2, others)
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


