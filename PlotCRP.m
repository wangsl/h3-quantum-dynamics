

clear all
%close all
clc

format long

H2eV = 27.21138505;

load('CRP.mat')

E = CRP.energies;
crp = CRP.CRP;


h_CRP = plot(E*H2eV, crp, 'b', 'LineWidth', 2.5, ...
	     'MarkerEdgeColor','r', ...
	     'MarkerFaceColor', 'y', 'MarkerSize', 3.5, ...
	     'YDataSource', 'crp');

grid on;

set(gca, 'xtick', [0.4:0.4:max(E)*H2eV]);
set(gca, 'ytick', [0.0:0.1:1.2]);

hold off

axis([min(E)*H2eV, max(E)*H2eV, -0.1, 1.1]);
