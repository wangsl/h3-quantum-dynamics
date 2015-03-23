
% Initial State Selected Cummulative Reaction Probabilities
% wave function in energy

function [ faiE, DfaiE, crp ] = ISSCRP2(r1, r2, theta, psi, ...
					time, CRP)

fprintf(' Calculate wave in energy\n');

dr2 = r2.dr;
n2DivSurf = CRP.n_dividing_surface;

nPoints = CRP.n_gradient_points;
n = int32((nPoints-1)/2);

fai = psi(1:2:end, (n2DivSurf-n):(n2DivSurf+n), :) + ...
      j*psi(2:2:end, (n2DivSurf-n):(n2DivSurf+n), :);

[ fai, Dfai ] = Gradient3(dr2, fai);

dt = time.time_step;
nt = time.steps;
t = double(nt)*dt;

E = CRP.energies; 

if nt == 1
  f = exp(j*E*t)*(dt/2);
else
  f = exp(j*E*t)*dt;
end

nE = numel(E);
n1 = size(fai, 1);
nTheta = size(fai, 3);

fai = reshape(fai, [n1*nTheta, 1]);
Dfai = reshape(Dfai, [n1*nTheta, 1]);

fai = fai*f;
Dfai = Dfai*f;

faiE = CRP.faiE + reshape(fai, [n1, nTheta, nE]);
DfaiE = CRP.DfaiE + reshape(Dfai, [n1, nTheta, nE]);

% Initial state selected reaction probability

eta2 = CRP.eta_sq;
dr1 = r1.dr;
mu2 = r2.mass;

% sum over r1
crp = sum(conj(faiE).*DfaiE);
crp = reshape(crp, [nTheta, nE]);

% sum over theta
crp = reshape(crp, [nTheta, nE]);
crp = theta.w'*crp;

crp = imag(crp)./eta2*(dr1/mu2);


