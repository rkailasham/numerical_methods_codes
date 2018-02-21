function [sol] = analytic(q)
%q refers to epsilon value
count=0;
f_actual = zeros(961,1);
u_actual = zeros(961,1);
h=1/30;
for i=1:31
    for j=1:31
        count = ((j-1)*31 + i);
        x_i = (i-1)*h;  %setting up the vector with actual solutions!!
                        %solving the given equation subject to the boundary conditions that
                        %u(x,y) = 0 at all the boundaries, we get u(x,y) =
                        %sin(pi*r*x)sin(pi*s*y).
        y_j = (j-1)*h;        
        f_actual(count)= q*sin(pi*x_i);% to generate the RHS of Au = f
        u_temp(i,j) = (sin(pi*1*x_i)*sin(pi*2*y_j));
        u_actual(count)= u_temp(i,j);
        %sol1(k,:) = Jac_Res(u_actual,f_actual,h,P(k),q,m,tol);        
                 
        end
end
var_count = 0;
lambda_o = (pi^2)*(2^2 + 1^2);

for lambda = lambda_o:0.1:51
    
    var_count = var_count + 1;
    
    [sol(var_count,:), J_u_inv, J_lambda] = newton_new(0.00001,30,lambda,0,u_actual,f_actual);
    delx_delp = J_u_inv * J_lambda;
    delx_delp = -delx_delp;
    u_actual = sol(var_count,:)' + (delx_delp)*(0.1);
    
end

end