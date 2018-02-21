%************************ INPUT PARAMETERS ***************************** %
%
% sol_branch0 : converged solution(say u0) at an initial lambda value.
% sol_branch1 : guess for a solution(say u1) at another lambda value.
% s           : step size in s.
% p0          : lambda value at sol_branch0
% p1          : initial gues for a lambda corresponding to sol_branch1
% q           : epsilon
% f_actual    : forcing function
% tol         : tolerance. Set at 10e-5
% no          : number of points on the unit grid in either x or y
%               direction. Set as 30 here.
%
%********************* OUTPUT VALUES ************************************%
%
% sol1        : converged value of the solution at given s.
% sol_lambda  : converged value of lambda at iven s.
% J_aug       : as explained on the first page of Appendix A.
% J_s         : as explained on the first page of Appendix A.
% P           : Permutation matrix at the end of the LU decomposition of
%               J_aug.
%
%************************************************************************%

function [sol1, sol_lambda, J_aug,J_s, P] = newton_ALC(sol_branch0,sol_branch1,s,p0,p1, q, f_actual,tol, no)
%p refers to lambda
%q refers to epsilon
h=1/no;
m = no + 1;
N = m^2;
count = 0;
order=[1:(N+1)];
u_actual0 = sol_branch0;
u_actual1 = sol_branch1;
norm_u=1;
norm_la = 1;
while(norm_u>tol && norm_la>tol) % to ensure convergence in both lambda AND u

[J_aug,R_aug, J_s] = Jac_Res_ALC(u_actual0,u_actual1,f_actual,h,p0,p1,q,m,s);
[L,U,P] = lu(J_aug);
sol = Substitute(lu(P*J_aug),order,(N+1),(P*R_aug));
del_u = sol(1:961);
del_la = sol(962);
norm_u = norm(del_u,2);
norm_la = norm(del_la);
u_t = u_actual1;
u_t = u_t + del_u;
u_actual1 = u_t;
p1 = p1 + del_la;
end
sol1 = u_actual1; 
sol_lambda = p1;
J_aug = lu(P*J_aug); %returns the LU decomposed Jacobian to the ALC method
end








    
    
        