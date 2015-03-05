
% $Id$

clear all
%close all

clc
format long

global masses
global UseLSTH

UseLSTH = true;

ElectronMass = 9.10938291E-31;
AMU = 1.66053892E-27;
MassAU = AMU/ElectronMass;

mH = 1.007825*MassAU;
mD = 2.01410178*MassAU;
mT = 3.0160492*MassAU;

masses = [ mH mH mH ];

R = linspace(0.2, 8.0, 128);
r = linspace(0.2, 8.0, 128);

n = 3;
m = 3;

k = 0;
for i = 1 : n
  for j = 1 : m
    k = k + 1;
    subplot(n, m, k)
    Theta = (k-1)*10/180.0*pi;
    
    tic
    V = H3PESJacobi(R, r, Theta);
    toc
    
    [ ~, hPES ] = contourf(R, r, V', [ -0.1:0.01:-0.001 0.001:0.01:0.3 ]);
    set(hPES, 'LineWidth', 0.75);
    set(hPES, 'LineColor', 'black');
    %colorbar('vert');
    colormap cool;
    axis square
  end
end

print -dpdf t1.pdf

return
