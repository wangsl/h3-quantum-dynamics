
% $Id$

function [ e, psi ] = H2WaveFunction(R, varargin)

n = R.n;
dr = R.dr;
mu = R.mass;
r = R.r;

if length(varargin) > 0 
  nVbs = [ varargin{:} ];
else
  nVbs = 0:1:n-1;
end

H = zeros(n);

for i = 1 : n
  H(i,i) = pi^2/3 - 1.0/(2*i*i);
  for j = 1 : i-1
    H(i,j) = (-1)^(i-j)*(2.0/(i-j)^2 - 2.0/(i+j)^2);
    H(j,i) = H(i,j);
  end
end

H = H/(2*mu*dr*dr);

V = H2PES(r);
V = reshape(V, [1, numel(V)]);

% Get diagonal elements
H(1:size(H,1)+1:end) = H(1:size(H,1)+1:end) + V;

[ vecs, energies ] = eig(H);

% BKMP2 PES
% vH2Min = -0.174495770896975;

nVbs = nVbs+1;
e = diag(energies);
e = e(nVbs);
psi = vecs(:, nVbs)/sqrt(dr);

return

