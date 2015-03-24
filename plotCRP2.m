
function [ ] = PlotCRP(CRP)

global H2eV

persistent has_CRP_plot 
persistent h_CRP 

E = CRP.energies;
crp = CRP.CRP;

if isempty(has_CRP_plot)
  
  has_CRP_plot = 1;
  
  figure(2)
  
  h_CRP = plot(E*H2eV, crp, 'b', 'LineWidth', 2.5, ...
               'MarkerEdgeColor','r', ...
               'MarkerFaceColor', 'y', 'MarkerSize', 3.5, ...
               'YDataSource', 'crp');
  
  grid on;
  
  set(gca, 'xtick', [0.4:0.1:max(E)*H2eV]);
  set(gca, 'ytick', [0.0:0.1:1.2]);
  
  hold off
  
  axis([min(E)*H2eV, max(E)*H2eV, -0.1, 1.1]);
end

refreshdata(h_CRP, 'caller');
drawnow
