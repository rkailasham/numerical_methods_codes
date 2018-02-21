function [sol1, J_u_inv,J_lambda] = newton_new(tol, no, p, q,sol_branch,f_actual)
%p refers to lambda
%q refers to epsilon
h=1/no;
syms x1;
syms x2;
m = no + 1;
N = m^2;
x_i=0;
y_j=0;
u_temp = zeros(m);
branch_count=0;
r=0;
s=0;
count = 0;

u_actual = sol_branch;

norm_u=1;
while(norm_u>tol)
[Jac_u,Resid_u, J_lambda] = Jac_Res_new(u_actual,f_actual,h,p,q,m,tol);
del_u_calc = Jac_u\Resid_u;
norm_u = norm(del_u_calc,2);
u_t = u_actual;
u_t = u_t + del_u_calc;
%disp(size(u_actual(k,:)));
%disp(size(del_u_calc));
u_actual = u_t;
end
sol1 = u_actual; 
J_u_inv = inv(Jac_u);
end








    
    
        