%********************** INPUT PARAMETERS *******************************%
% n  - Number of points at which the solution needs to be evaluated     %
% x0 - Initial condition for x                                          %
% y0 - Initial condition for y                                          %
% tspan  - Timespan over which solutions need to be evaluated           %
% a,b,c,d,p,q - Constants as defined in the problem statement.          %
%                                                                       %
%********************** OUTPUT VALUES **********************************%
%                                                                       %
% solx - The solution x(t) calculated using implicit Euler method       %
% soly - The solution y(t) calculated using implicit Euler method       %
% time - Points on the time axis at which solutions have been           %
%        calculated                                                     %
%                                                                       %
%************************************************************************


function [solx,soly,t] = ImpODE(n,x0,y0,tspan,a,b,c,d,p,q)
syms x;
syms y;
dx_dt = (a-b*y)*x - p*x^2;
dy_dt = (c*x-d)*y - q*y^2;

solx(1) = x0;
soly(1) = y0;

h=tspan/n;
t(1) = 0;
tol = (h^2); % Tolerance for Newton's iteration set as h^2.
for i=1:n
      
    %Generating initial guesses using Predictor-corrector
    fx = func_eval(dx_dt,solx(i),soly(i));
    fy = func_eval(dy_dt,solx(i),soly(i));
    
    x_init = solx(i) + h*fx;
    y_init = soly(i) + h*fy;
    
    t(i+1) = i*h;
    [solx(i+1),soly(i+1)] = newton(tol,x_init,y_init,solx(i),soly(i),h,a,b,c,d,p,q);
    
end
end
