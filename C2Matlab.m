
% $Id$

function [] = C2Matlab(r1, r2, theta, pot, psi, time, options, others)

fprintf(' From C2Matlab\n')

ISSCRP(r1, r2, theta, psi, others.CRP, others.divSurf, time);

if mod(time.steps, 20) == 0
  PlotPotWave(r1, r2, pot, psi)
end
