
% $Id$

function [ V ] = H3PESJacobi(R, r, Theta)

global masses

f2 = masses(2)/(masses(2)+masses(3));
f3 = 1 - f2;

nR = length(R);
nr = length(r);
nTheta = length(Theta);

VMax = 4.0;
rMin = 0.4;

V = zeros(nR, nr, nTheta);
V = V';
V(:,:,:) = VMax;

[x, y, z] = meshgrid(R, r, Theta);

r23 = y;
r12 = sqrt((f2*y).^2 + x.^2 - 2*f2*y.*x.*cos(z));
r13 = sqrt((f3*y).^2 + x.^2 + 2*f3*y.*x.*cos(z));

Grids = find(r12>rMin & r13>rMin & r23>rMin);
V(Grids) = H3PES(r12(Grids), r23(Grids), r13(Grids));

V = V';

return


