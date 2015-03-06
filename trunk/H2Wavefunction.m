
% $Id$

function [ e, psi ] = H2Wavefunction(R, nVb)

n = R.n;
dr = R.dr;
mu = R.mass;;
r = R.r;

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

% Get diagonal elements
H(1:size(H,1)+1:end) = H(1:size(H,1)+1:end) + V';

[ vecs, energies ] = eig(H);

%vH2Min = -0.174495770896975;
%e = diag(energies) + vH2Min;
%psi = vecs/sqrt(dr);
%return

e = energies(nVb+1, nVb+1);
psi = transpose(vecs(:, nVb+1)/sqrt(dr));

return

