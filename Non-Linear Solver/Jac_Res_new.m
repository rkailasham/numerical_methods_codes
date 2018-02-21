function[J_u,R_u,J_lambda] = Jac_Res_new(u_input, f_input, h, p,q,m, tol)
%p refers to lambda
%q refers to epsilon
%m = (n+1) i.e number of points for fin diff grid +1
u_input = u_input';
f_input = f_input';
N = m^2;
J_u = zeros(N);
J_lambda = zeros(N,1);
norm_u = 1;
del_u_calc = zeros(N,1);
R_u = zeros(N,1);

count=0;


  
for i=1:m
    for j=1:m
        count = ((j-1)*m + i);
        
        % setting up co-efficients for the nodes at the edges
         if(i==1 && j==1) %point (1,1)
            J_u(count,count) = 1;
            J_lambda(count) = 0;
            R_u(count) = 0;
            
        elseif (i==1 && j==m) %point(1,m)
           J_u(count,count) = 1;
           J_lambda(count) = 0;
           R_u(count) = 0;
        elseif (i==m && j==1) %point(m,1)
            J_u(count,count) = 1;
            J_lambda(count) = 0;
            R_u(count) = 0;
         elseif (i==m && j==m) %point(m,m)
            J_u(count,count) = 1;
            J_lambda(count) = 0;
            R_u(count) = 0;
        % end of setting up co-efficients for nodes
        
        %setting up coefficients for points ON the boundary line
        %using central difference for the convection term (first
        %derivative) and appropriate method for second derivative.
        
        elseif (i==1 && i~=j && j~=m)
           J_u(count,count) = 1;
           J_lambda(count) = 0;
           R_u(count) = 0; 
        elseif (i==m && i~=j && j~=1)
            J_u(count,count) = 1;
            J_lambda(count) = 0;
            R_u(count) = 0;
        elseif(j==1 && i~=j && i~=m)
            J_u(count,count) = 1;
            J_lambda(count) = 0;
            R_u(count) = 0;
        elseif(j==m && i~=j && i~=1)
           J_u(count,count) = 1;
           J_lambda(count) = 0;
           R_u(count) = 0;
         else
            R_u(count) = u_input(count)*((-4/h^2)+ p) + p*u_input(count)^2 + (1/h^2) * (u_input(count-1) + u_input(count+1) + u_input(count+m) + u_input(count-m)) - f_input(count); 
            J_u(count,count) = -4/(h^2) + p + 2*p*u_input(count);
            J_u(count,(count-1)) = 1/(h^2);
            J_u(count,(count+1)) = 1/(h^2);
            J_u(count,(count-m)) = 1/(h^2);
            J_u(count,(count+m)) = 1/(h^2);
            
            J_lambda(count) = u_input(count) + (u_input(count))^2;
        end
     end
end
R_u = -R_u;

end

