
% $Id$

function [] = C2Matlab() %r1, r2, theta, pot, psi, time, options, others)

global H3Data

fprintf(' From C2Matlab\n')

ISSCRP() %r1, r2, theta, psi, others.CRP, others.divSurf, time);

if mod(H3Data.time.steps, 20) == 0
  PlotPotWave(H3Data.r1, H3Data.r2, H3Data.pot, H3Data.psi)
end
