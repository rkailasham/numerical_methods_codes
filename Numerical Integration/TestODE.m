
%********************** INPUT PARAMETERS *******************************%
% n  - Number of points at which the solution needs to be evaluated     %
% x0 - Initial condition for x                                          %
% y0 - Initial condition for y                                          %
% T  - Timespan over which solutions need to be evaluated               %
% a,b,c,d,p,q - Constants as defined in the problem statement.          %
%                                                                       %
%********************** OUTPUT VALUES **********************************%
%                                                                       %
% normx - Absolute value of the difference between x(at a point)        %
%         calculated using the two methods, normalised by the           %
%         number of points                                              %
% normy - Absolute value of the difference between y(at a point)        %
%         calculated using the two methods, normalised by the           %
%         number of points                                              %
%************************************************************************

function[normx,normy] = TestODE(n,x0,y0,T,a,b,c,d,p,q)

solx = zeros(n+1,1);
soly = zeros(n+1,1);
time = zeros(n+1,1);

solx(1) = x0;
soly(1) = y0;
time(1) = 0;
h=T/n;

fx(1) = (a-b*soly(1))*solx(1) - p*(solx(1)^2);
fy(1) = (c*solx(1)-d)*soly(1) - q*(soly(1)^2);

%using Euler's method to calculate the solution at the first step
solx(2) = solx(1) + h*fx(1);
soly(2) = soly(1) + h*fy(1);

fx(2) = (a-b*soly(2))*solx(2) - p*(solx(2)^2);
fy(2) = (c*solx(2)-d)*soly(2) - q*(soly(2)^2);


solx(3) = solx(2) + (h/2)*(3*fx(2)-fx(1));
soly(3) = soly(2) + (h/2)*(3*fy(2)-fy(1));

for i=3:n
    fx(i) = (a-b*soly(i))*solx(i) - p*(solx(i)^2);
    fy(i) = (c*solx(i)-d)*soly(i) - q*(soly(i)^2);
    time(i+1) = i*h;
    %using Explicit 2nd order Adams-Bashforth for error analysis.
    solx(i+1) = solx(i) + (h/2)*(3*fx(i)-fx(i-1));
    soly(i+1) = soly(i) + (h/2)*(3*fy(i)-fy(i-1));
end

[solxE,solyE,time] = ImpODE(n,x0,y0,T,a,b,c,d,p,q); 

%this line can be changed to ExpODE to calculate error with respect to
%Explicit Euler.

%evaluating errors at the midway point along the time axis.
normx = norm (solx(ceil(n/2)) -solxE(ceil(n/2)))/n;
normy = norm (soly(ceil(n/2)) -solyE(ceil(n/2)))/n;

end