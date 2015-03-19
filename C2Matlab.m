
% $Id$

function [] = C2Matlab(r1, r2, theta, pot, psi, time, options)

fprintf(1, ' test from C2Matlab\n')

if mod(time.steps, 10) == 0
  PlotPotWave(r1, r2, pot, psi)
end
