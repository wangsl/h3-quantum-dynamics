
% $Id$

function [ v, g ] = Gradient3(dx, f)

n = size(f, 2);

if n == 3 
  [ v, g ] = Gradient_3(dx, f);
elseif n == 5
  [ v, g ] = Gradient_5(dx, f);
elseif n == 7
  [ v, g ] = Gradient_7(dx, f);
elseif n == 9
  [ v, g ] = Gradient_9(dx, f);
elseif n == 11 
  [ v, g ] = Gradient_11(dx, f);
else
  error('Unsupport input arry size')
end

return

function [ v, g ] = Gradient_3(dx, f) 
v = f(:,2,:);
g = f(:,3,:)-f(:,1,:);
g = g/(2*dx);
return

function [ v, g ] = Gradient_5(dx, f) 
v = f(:,3,:);
g = -(f(:,5,:)-f(:,1,:)) + 8*(f(:,4,:)-f(:,2,:));
g = g/(12*dx);
return

function [ v, g ] = Gradient_7(dx, f) 
v = f(:,4,:);
g = (f(:,7,:)-f(:,1,:)) - 9*(f(:,6,:)-f(:,2,:)) + 45*(f(:,5,:)-f(:,3,:));
g = g/(60*dx);
return

function [ v, g ] = Gradient_9(dx, f)
v = f(:,5,:);
g = -3*(f(:,9,:)-f(:,1,:)) + 32*(f(:,8,:)-f(:,2,:)) ...
    - 168*(f(:,7,:)-f(:,3,:)) + 672*(f(:,6,:)-f(:,4,:));
g = g/(840*dx);
return

function [ v, g ] = Gradient_11(dx, f)
v = f(:,6,:);
g = 2*(f(:,11,:)-f(:,1,:)) - 25*(f(:,10,:)-f(:,2,:)) ...
    + 150*(f(:,9,:)-f(:,3,:)) - 600*(f(:,8,:)-f(:,4,:)) ...
    + 2100*(f(:,7,:)-f(:,5,:));
g = g/(2520*dx);
return

