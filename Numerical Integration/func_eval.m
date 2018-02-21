function [sol] = func_eval(expression,a,b)
syms x;
syms y;
sol = eval(subs(expression,{x,y},{a,b}));
end
