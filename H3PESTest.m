
% $Id$

clear all
%close all

clc
format long

global masses
global UseLSTH

% UseLSTH = true;

ElectronMass = 9.10938291E-31;
AMU = 1.66053892E-27;
MassAU = AMU/ElectronMass;

mH = 1.007825*MassAU;
mD = 2.01410178*MassAU;
mT = 3.0160492*MassAU;

masses = [ mH mH mH ];

R = linspace(0.4, 9.0, 512);
r = linspace(0.4, 7.0, 512);
Theta = linspace(0.0, pi, 19);

tic
V = H3PESJacobi(R, r, Theta, masses);
toc

n = 3;
m = 3;

k = 0;
for i = 1 : n
  for j = 1 : m
    k = k + 1;
    subplot(n, m, k)
    
    [ ~, hPES ] = contourf(R, r, V(:,:,k)', [ -0.1:0.01:-0.001 ...
		    0.001:0.01:0.20 ]);
    set(gca, 'xtick', [1:2:max(R)]);
    set(gca, 'ytick', [1:2:max(r)]);
    set(hPES, 'LineWidth', 0.75);
    set(hPES, 'LineColor', 'black');
    %colorbar('vert');
    colormap cool;
    %axis square
  end
end

print -dpdf H3PESJacobi.pdf

return
