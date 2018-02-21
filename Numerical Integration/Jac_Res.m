%************************** INPUT PARAMETERS **************************%
%                                                                      %
% x_init - Initial guess for x at the given point                      %
% y_init - Initial guess for y at the given point                      %
% x_prev - converged value of x at the previous point                  %
% y_prev - converged value of y at the previous point                  %
% h      - step size in time                                           %
% a,b,c,d,p,q - Constants as defined in the problem statement          %
%                                                                      %
%************************** OUTPUT VALUES *****************************%
%                                                                      %
% J -    Jacobian matrix calculated at (x_init,y_init)                 %
% R -    Residual vector calculated at (x_init,y_init)                 %
%                                                                      %
%**********************************************************************%

function[J,R] = Jac_Res(x_init,y_init,x_prev,y_prev,h,a,b,c,d,p,q)

R1 = x_init - x_prev - h*((a-b*y_init)*x_init - p*x_init^2); %eqn for dx_dt
R2 = y_init - y_prev -h*((c*x_init-d)*y_init - q*y_init^2);  %eqn for dy_dt

J1x = 1 - h* ( (a-b*y_init) - 2*p*x_init);
J1y = -h*(-b*x_init);

J2x = -h*(c*y_init);
J2y = 1 - h*( (c*x_init-d) - 2*q*y_init);

J(1,1) = J1x;
J(1,2) = J1y;
J(2,1) = J2x;
J(2,2) = J2y;

R = -1*cat(1,R1,R2);

end

