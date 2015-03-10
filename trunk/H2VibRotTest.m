
% $Id$

clear all
close all
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

r.n = int64(4096);
r.r = linspace(0.4, 16.0, r.n);
r.dr = r.r(2) - r.r(1);
r.mass = mH/2;

R = r.r;
mu = r.mass;

%jRot = 0;
%plot(R, H2PES(R)' + jRot*(jRot+1)./(2*mu*R.^2));
%axis([min(R) max(R) 0 0.2]);
%return

energy_levels = 0:1:15; %[ 0 1 6 13 14 15 16 ];

[ e, psi ] = H2VibRotWaveFunction(r, 5, energy_levels);

plot(r.r, psi, 'LineWidth', 2);
legend(arrayfun(@num2str, energy_levels, 'UniformOutput', false));

print -dpdf H2-VibRot.pdf

return

  -0.156641476303803
  -0.138080296682007
  -0.120583361684047
  -0.104130771802557
  -0.088711613237419
  -0.074324899327577
  -0.060980960361705
  -0.048703468324473
  -0.037532360745080
  -0.027528054100074
  -0.018777594858582
  -0.011403986746862
  -0.005581378529287
  -0.001563620937494
   0.000166134772494
   0.000246533536189
