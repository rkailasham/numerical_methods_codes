%***********************************************************************%
% This function takes in 2 solutions and computes the error of the      %
% coarser solution (the one calculated using a smaller number of        %
% elements) with respect to the finer solution. The error is normalised %
% by the number of points in the coarser solution.                      %
%                                                                       %
%***********************************************************************%

function [err] = CalcError(sol_fine,sol_coarse)
l1 = length(sol_fine);
l2 = length(sol_coarse);
t= zeros(l2,1);
fact = ((l1-1)/(l2-1));
disp(fact);
for i=1:l2
    pos = 1 + (i-1)*fact;
    t(i) = ((sol_fine(pos) - sol_coarse(i)));
end
err = norm(t)/(l2);
end