%******************* INPUT PARAMETERS *********************************** %
% q    : value of epsilon at which analytic continuation is to be carried % 
%        out                                                              %
% mode : Will perform analytic continuation in epsilon if this parameter  %
%        is entered as 'epsilon'.                                         %
%                                                                         %
%****************** OUTPUT VALUES ****************************************%
%
% sol_la      :  Returns the family of solutions that have been stepped up 
%                in lambda at initial value of epsilon (usually 0).       %
% sol_ep_temp : Contains the family of solutions at a given epsilon value %
%               for one particular lambda value.                          %
% sol_ep      : Contains the family of solution over the entire range of  %
%                lambda for eosilon value not equal to zero.              %
% norm_la     : contains the 2-norm of each solution in sol_la            %
% norm_ep     : contains the 2-norm of each vector in sol_ep              % 
%                                                                         %
%*************************************************************************%

function [sol_la, sol_ep_temp, sol_ep, norm_la,norm_ep] = analytic(q,mode)
%q refers to epsilon value
count=0;
sol_ep_temp=zeros(11,961);
sol_ep = zeros(11,961);
norm_ep = zeros(1,10);
norm_la = zeros(1,10);
f_actual = zeros(961,1);
u_actual = zeros(961,1);
h=1/30;
for i=1:31
    for j=1:31
        count = ((j-1)*31 + i);
        x_i = (i-1)*h;  %setting up the vector with actual solutions!!
                        %solving the given equation subject to the boundary conditions that
                        %u(x,y) = 0 at all the boundaries, we get u(x,y) =
                        %sin(pi*m*x)sin(pi*n*y).
        y_j = (j-1)*h;        
        f_actual(count)= sin(pi*x_i);% to generate the RHS of Au = f
        u_temp(i,j) =1 * (sin(pi*1*x_i)*sin(pi*1*y_j));
        u_actual(count)= u_temp(i,j);
    end
end
  
var_count = 0;
lambda_o = (pi^2)*(1^2 + 1^2) -0.15;
order=[1:961];
u_initial = u_actual;

for lambda = lambda_o:-0.1:10
    
    var_count = var_count + 1;
    
    [sol_la(var_count,:), J_P, P, J_lambda,J_epsilon] = newton_new(0.00001,30,lambda,q,u_actual,f_actual);
    delx_delp = Substitute((J_P),order,961,(P*J_lambda));
    u_actual = sol_la(var_count,:)' + (delx_delp)*(-0.1);
    
end
u_actual = sol_la(var_count,:)';
lambda_o = 10.0892;
var_count=0;

for lambda = lambda_o:0.1:60
    
    var_count = var_count + 1;
    
    [sol_la(var_count,:), J_P, P, J_lambda,J_epsilon] = newton_new(0.00001,30,lambda,q,u_actual,f_actual);
    delx_delp = Substitute((J_P),order,961,(P*J_lambda));
    norm_la(var_count) = 1*norm(sol_la(var_count,:));
   
    u_actual = sol_la(var_count,:)' + (delx_delp)*(0.1);
    
end

if(strcmp(mode,'epsilon'))
%implementing analytic continuation in epsilon    
var_count=0;
u_actual = sol_la(1,:)';


for epsilon =q:0.1:1
    var_count = var_count+1;
    
    [sol_ep_temp(var_count,:), J_P, P, J_lambda,J_epsilon] = newton_new(0.00001,30,lambda_o,epsilon,u_actual,f_actual);
    delx_delq = Substitute((J_P),order,961,(P*J_epsilon));
    %disp(size(delx_delp));
    %disp(size(sol(var_count,:)'));
    u_actual = sol_ep_temp(var_count,:)' + (delx_delq)*(0.1);
end


u_actual = sol_ep_temp(11,:)';
q=1.0;
var_count=0;
lambda_o = 10.0892;

for lambda = lambda_o:0.1:60
    
    var_count = var_count + 1;
    
    [sol_ep(var_count,:), J_P, P, J_lambda,J_epsilon] = newton_new(0.00001,30,lambda,q,u_actual,f_actual);
    delx_delp = Substitute((J_P),order,961,(P*J_lambda));
    norm_ep(var_count) = norm(sol_ep(var_count,:));
    u_actual = sol_ep(var_count,:)' + (delx_delp)*(0.1);
    
end
end
end