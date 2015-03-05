
% $Id$

function [ V ] = H2PES(r)

r1 = reshape(r, [numel(r), 1]);
r2 = zeros(size(r1));
r2(:) = 100.0;

V = H3PES(r1, r2, r1+r2);
return
