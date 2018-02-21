%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function accepts f, afunction of x and y, and the x and y values at%
% which the value of the function has to be evaluated.                    %
% The function returns the value of f at the pecified values.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol] = func_eval(expression,a,b)
syms x;
syms y;
sol = eval(subs(expression,{x,y},{a,b}));
end
