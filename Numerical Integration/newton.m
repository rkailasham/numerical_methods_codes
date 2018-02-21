%************************** INPUT PARAMETERS **************************%
%                                                                      %
% tol    - Tolerance value. Depends on the step-size,h                 %
% x_init - Initial guess for x at the given point                      %
% y_init - Initial guess for y at the given point                      %
% x_prev - converged value of x at the previous point                  %
% y_prev - converged value of y at the previous point                  %
% h      - step size in time                                           %
% a,b,c,d,p,q - Constants as defined in the problem statement          %
%                                                                      %
%************************** OUTPUT VALUES *****************************%
%                                                                      %
% solx - converged value of x at the given point                       %
% soly - converged value of y at the given point                       %
%                                                                      %
%**********************************************************************%

function [solx,soly] = newton(tol,x_init,y_init,x_prev,y_prev,h,a,b,c,d,p,q)

norm_x=1;
norm_y = 1;

while(norm_x>tol && norm_y > tol)
[Jac_u,Resid_u] = Jac_Res(x_init,y_init,x_prev,y_prev,h, a,b,c,d,p,q);
del_u_calc = Jac_u\Resid_u;
del_x = (del_u_calc(1));
del_y = (del_u_calc(2));
norm_x = norm(del_x,2);
norm_y = norm(del_y,2);
x_init = x_init + del_x;
y_init = y_init + del_y;
end

solx = x_init;  %this x_init contains the converged value of x at point (n+1)
soly = y_init;  %this y_init contains the converged value of y at point (n+1)

end








    
    
        