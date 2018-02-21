%************************ INPUT PARAMETERS ******************************%
%
% u_input0 : one value of solution.
% u_input1 : another value of the solution.
% f_input  : forcing function
% p0       : lambda value corresponding to u_input0
% p1       : lambda value corresponding to u_input1
% q        : epsilon value at which ALC is performed
% m        : (n+1) i.e number of points for finite difference grid +1
% s        : step size in s. 
%
%************************ OUTPUT VALUES *********************************%
%
% J_aug : as explained on the first page of Appendix A.
% R_aug : as explained on the first page of Appendix A.
% J_s   : as explained on the first page of Appendix A.
%
%************************************************************************%


function[J_aug,R_aug,J_s] = Jac_Res_ALC(u_input0,u_input1, f_input, h, p0,p1,q,m,s)

u_input0 = u_input0';
u_input1 = u_input1';
f_input = f_input';
N = m^2;
J_u = zeros(N);
J_lambda = zeros(N,1);
J_epsilon = zeros(N,1);
J_s = zeros(N+1,1);
R_u = zeros(N,1);
eta_u = zeros(1,N);


count=0;

for i=1:m
    for j=1:m
        count = ((j-1)*m + i);
        
        % setting up co-efficients for the nodes at the edges
         if(i==1 && j==1) %point (1,1)
            J_u(count,count) = 1;
            
            R_u(count) = 0;
            
        elseif (i==1 && j==m) %point(1,m)
           J_u(count,count) = 1;
          
           R_u(count) = 0;
        elseif (i==m && j==1) %point(m,1)
            J_u(count,count) = 1;
            
            R_u(count) = 0;
         elseif (i==m && j==m) %point(m,m)
            J_u(count,count) = 1;
           
            R_u(count) = 0;
        % end of setting up co-efficients for nodes
        
        %setting up coefficients for points ON the boundary line
        %using central difference.
        
        elseif (i==1 && i~=j && j~=m)
           J_u(count,count) = 1;
          
           R_u(count) = 0; 
        elseif (i==m && i~=j && j~=1)
            J_u(count,count) = 1;
            
            R_u(count) = 0;
        elseif(j==1 && i~=j && i~=m)
            J_u(count,count) = 1;
            
            R_u(count) = 0;
        elseif(j==m && i~=j && i~=1)
           J_u(count,count) = 1;
           R_u(count) = 0;
           
         else
            R_u(count) = u_input1(count)*((-4/h^2)+ p1) + p1*u_input1(count)^2 + (1/h^2) * (u_input1(count-1) + u_input1(count+1) + u_input1(count+m) + u_input1(count-m)) - q*f_input(count); 
            J_u(count,count) = -4/(h^2) + p1 + 2*p1*u_input1(count);
            J_u(count,(count-1)) = 1/(h^2);
            J_u(count,(count+1)) = 1/(h^2);
            J_u(count,(count-m)) = 1/(h^2);
            J_u(count,(count+m)) = 1/(h^2);
            
            J_lambda(count) = u_input1(count) + (u_input1(count))^2;
            eta_u(count) = -2*(u_input1(count) - u_input0(count));
        end
     end
end

eta = (s)^2 - (norm(u_input1-u_input0))^2 - (p1-p0)^2;
eta_lambda = -2*(p1-p0);
eta_u = cat(2,eta_u,eta_lambda);
J_u = cat(2,J_u,J_lambda);
J_aug = cat(1,J_u,eta_u);
R_aug = -1* cat(1,R_u,eta);

%%J_aug and J_s are needed for the arc-length continuation method
J_s(N+1) = 2*s;
J_s = -1* J_s;


end

