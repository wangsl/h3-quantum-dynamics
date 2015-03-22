
% $Id$

function [] = C2Matlab() 

global H3Data

fprintf(' From C2Matlab\n')

[ H3Data.CRP.faiE, H3Data.CRP.DfaiE, H3Data.CRP.CRP ] = ...
    ISSCRP2(H3Data.r1, H3Data.r2, H3Data.theta, ...
	    H3Data.psi, H3Data.divSurf, ...
	    H3Data.time, H3Data.CRP);

if mod(H3Data.time.steps, 20) == 0
  PlotPotWave(H3Data.r1, H3Data.r2, H3Data.pot, H3Data.psi)
  PlotCRP(H3Data.CRP);
end
