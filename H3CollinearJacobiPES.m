
% $Id$

function [ V ] = H3CollinearJacobiPES(R, r) 

global masses

nR = length(R);
nr = length(r);

f = masses(3)/(masses(2)+masses(3));

V = zeros(nR, nr);

VMax = 4.0;
rMin = 0.4;

V(:,:) = VMax;

for j = 1 : nr
  r2 = r(j);
  if r2 < rMin; continue; end
  
  r1 = R - f*r2;
  r1Index = find(r1 >= rMin);
  
  r1Tmp = r1(r1Index);
  r2Tmp = zeros(size(r1Tmp));
  r2Tmp(:) = r2;
  
  V(r1Index, j) = BKMP2(r1Tmp, r2Tmp, r1Tmp+r2Tmp);
end

return

% $$$ 
% $$$ V = zeros(size(V));
% $$$ 
% $$$ for i = 1 : nR
% $$$   for j = 1 : nr
% $$$     r2 = r(j);
% $$$     r1 = R(i) - f*r2;
% $$$     if r2 < rMin || r1 < rMin
% $$$       V(i,j) = VMax;
% $$$     else
% $$$       V(i,j) = BKMP2(r1, r2, r1+r2);
% $$$     end
% $$$   end
% $$$ end
% $$$ 
% $$$ 
