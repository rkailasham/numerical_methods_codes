%********************** INPUT PARAMETERS *******************************%
% n  - Number of points at which the solution needs to be evaluated     %
% x0 - Initial condition for x                                          %
% y0 - Initial condition for y                                          %
% T  - Timespan over which solutions need to be evaluated               %
%  a,b,c,d,p,q - Constants as defined in the problem statement.         %
%                                                                       %
%********************** OUTPUT VALUES **********************************%
%
% solx - The solution x(t) calculated using explicit Euler method       %
% soly - The solution y(t) calculated using explicit Euler method       %
% time - Points on the time axis at which solutions have been           %
%        calculated                                                     %
%                                                                       %
%************************************************************************

function[solx,soly,time] = ExpODE(n,x0,y0,T,a,b,c,d,p,q)

solx(1) = x0;
soly(1) = y0;
time(1) = 0;
h=T/n;
for i=1:n
    fx(i) = (a-b*soly(i))*solx(i) - p*(solx(i)^2);
    fy(i) = (c*solx(i)-d)*soly(i) - q*(soly(i)^2);
    time(i+1) = i*h;
    solx(i+1) = solx(i) + h*fx(i);
    soly(i+1) = soly(i) + h*fy(i);
end
end
