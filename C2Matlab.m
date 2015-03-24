
% $Id$

function [] = C2Matlab() 

global H3Data

fprintf(' From C2Matlab\n')

return

if mod(H3Data.time.steps, 200) == 0
  PlotPotWave(H3Data.r1, H3Data.r2, H3Data.pot, H3Data.psi)
end

if mod(H3Data.time.steps, 100) == 0
  PlotCRP(H3Data.CRP);
end

return

[ H3Data.CRP.faiE, H3Data.CRP.DfaiE, H3Data.CRP.CRP ] = ...
    ISSCRP2(H3Data.r1, H3Data.r2, H3Data.theta, ...
	    H3Data.psi, ...
	    H3Data.time, H3Data.CRP);

if mod(H3Data.time.steps, 20) == 0
  PlotCRP(H3Data.CRP);
end
