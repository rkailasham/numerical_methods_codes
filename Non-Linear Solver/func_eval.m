function [sol] = func_eval(expression,a,b)
syms x1;
syms x2;
sol = eval(subs(expression,{x1,x2},{a,b}));
end
