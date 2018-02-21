%************************** INPUT PARAMETERS ****************************%
%
% tol        : Tolerance value. Set as 10e-5
% no         : No. of points on the unit grid along x or y direction.
%              Set as 30.
% p          : value of lambda
% q          : value of epsilon
% sol_branch : initial guess
%f_actual    : forcing function
%
%************************* OUTPUT VALUES ********************************%
%
% sol1      : converged solution
% J_P       : LU decomposed Jacobian
% P         : Permutation matrix obtained during LU decomposition.
% J_lambda  : as explained on the first page of Appendix A.
% J_epsilon : as explained on the first page of Appendix A.
%
%************************************************************************%


function [sol1, J_P,P,J_lambda,J_epsilon] = newton_new(tol, no, p, q,sol_branch,f_actual)
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
order=[1:N];
u_actual = sol_branch;

norm_u=1;
while(norm_u>tol)
[Jac_u,Resid_u, J_lambda,J_epsilon] = Jac_Res_new(u_actual,f_actual,h,p,q,m);
[L,U,P] = lu(Jac_u);
del_u_calc = Substitute(lu(P*Jac_u),order,N,(P*Resid_u));
norm_u = norm(del_u_calc,2);
u_t = u_actual;
u_t = u_t + del_u_calc;
u_actual = u_t;
end
sol1 = u_actual;   %returns the converged solution to the analytic routine
J_P = lu(P*Jac_u); %returns the LU decomposed Jacobian to the analytic routine
end








    
    
        