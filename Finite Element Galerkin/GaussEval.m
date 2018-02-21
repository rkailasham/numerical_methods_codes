%***********************************************************************%
% This function uses the 3-point Gaussian quadrature to numerically     %
% evaluate the given function between the intervals -1 and +1           %
%***********************************************************************%

function [val] = GaussEval(expression)
syms e;

gauss_weights = [(5/9.),(8/9.),(5/9.)];
gauss_points  = [-1*sqrt(3/5),0,sqrt(3/5)];

a = 0;
for i=1:3
    a = a+ (gauss_weights(i)* eval(subs(expression,e,gauss_points(i))));
end
val = a;
end

