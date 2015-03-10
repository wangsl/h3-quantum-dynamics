
% $Id$

function [] = H2WaveFunctionAnalysis(r, phi, mu, J)

n = length(r);
dr = r(2) - r(1);

L = n*dr;
k = (2*pi/L)*[0:n/2 (-n/2+1):(-1)];

% avergaed potential energy
V = H2PES(r);
aV = sum(phi.^2.*V)*dr

% averaged rotational energy
V = J*(J+1)./(2*mu*r.^2);
V = V';
aRot = sum(phi.^2.*V)*dr

% averaged kinetic energy
f = fft(phi);
aT = sum(conj(f).*k'.^2.*f)*dr/n/(2*mu)

% total H2 energy
aT + aRot + aV

