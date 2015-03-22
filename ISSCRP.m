
% Initial State Selected Cummulative Reaction Probabilities
% (ISSCRP)

function [ ] = ISSCRP()
%r1, r2, theta, psi, CRPIn, divSurf, ...
%time)

global H3Data
global H2eV

fprintf(' Calculate wave in energy\n');

persistent has_CRP_plot 
persistent h_CRP 

%persistent CRP

%if isempty(CRP)
%  CRP = CRPIn;
%end

dr2 = H3Data.r2.dr;

n2DivSurf = H3Data.divSurf.n;

nPoints = H3Data.CRP.n_gradient_points;
n = int32((nPoints-1)/2);

fai = H3Data.psi(1:2:end, (n2DivSurf-n):(n2DivSurf+n), :) + ...
      j*H3Data.psi(2:2:end, (n2DivSurf-n):(n2DivSurf+n), :);

[ fai, Dfai ] = Gradient3(dr2, fai);

dt = H3Data.time.time_step;
nt = H3Data.time.steps;
t = double(nt)*dt;

E = H3Data.CRP.energy;
eH2 = H3Data.CRP.eH2;

if nt == 1
  f = exp(j*(E+eH2)*t)*(dt/2);
else
  f = exp(j*(E+eH2)*t)*dt;
end

nE = numel(E);
n1 = size(fai, 1);
nTheta = size(fai, 3);

fai = reshape(fai, [n1*nTheta, 1]);
Dfai = reshape(Dfai, [n1*nTheta, 1]);

fai = fai*f;
Dfai = Dfai*f;

H3Data.CRP.faiE = H3Data.CRP.faiE + reshape(fai, [n1, nTheta, nE]);
H3Data.CRP.DfaiE = H3Data.CRP.DfaiE + reshape(Dfai, [n1, nTheta, nE]);

if mod(nt, 20) ~= 0
  return
end

faiE = H3Data.CRP.faiE;
DfaiE = H3Data.CRP.DfaiE;
eta2 = H3Data.CRP.etaSq;
dr1 = H3Data.r1.dr;
mu2 = H3Data.r2.mass;

% sum over r1
crp = sum(conj(faiE).*DfaiE);
crp = reshape(crp, [nTheta, nE]);

% sum over theta
crp = reshape(crp, [nTheta, nE]);
crp = H3Data.theta.w'*crp;

H3Data.CRP.CRP = imag(crp)./eta2*(dr1/mu2);

crp = H3Data.CRP.CRP;

if isempty(has_CRP_plot)
  
  has_CRP_plot = 1;
  
  figure(2)
  
  h_CRP = plot((E+eH2)*H2eV, crp, 'b', 'LineWidth', 2.5, ...
               'MarkerEdgeColor','r', ...
               'MarkerFaceColor', 'y', 'MarkerSize', 3.5, ...
               'YDataSource', 'crp');
  
  grid on;
  
  set(gca, 'xtick', [0.4:0.1:max(E+eH2)*H2eV]);
  set(gca, 'ytick', [0.0:0.1:1.2]);
  
  hold off
  
  axis([min(E+eH2)*H2eV, max(E+eH2)*H2eV, -0.1, 1.1]);
end

refreshdata(h_CRP, 'caller');
drawnow

