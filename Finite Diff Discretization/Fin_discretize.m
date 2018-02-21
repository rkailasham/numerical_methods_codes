%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function was used to solve the no-convection case and to generate %
% reference grids for the convection case.                               %
%                                                                        %
%%%%%%%%%%%%%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
% n      : Number of grid points to be used                              %
% val_x  : The x-component of velocity.                                  %
% val_y  : The y-component of velocity.                                  %
% f      : Forcing function                                              %
%                                                                        %
%%%%%%%%%%%%%%%% FUNCTION OUTPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
% sol1   : Contains the solution u of the system Au = f                  %
% sol2   : Contains the matrix version of the solution u                 %
% sol3   : Contains the solution error nomalised by the number of grid   %
%          points (i.e (n+1) )                                           %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [sol1, sol2, sol3] = Fin_discretize(n,val_x,val_y,f)
h = 1/n;
N = (n+1)^2;
m = (n+1);
A = zeros(N,N);
x_i=0;
y_j=0;
syms x;
syms y;
vx=0;
vy=0;
u_temp = zeros(n+1);
u_actual = zeros(N,1);
f_actual = zeros(N,1);
count=0;
for i=1:m
    for j=1:m
        count = ((j-1)*(n+1) + i); %Converting the double-index system to a single index
        x_i = (i-1)*h;  
        y_j = (j-1)*h;
        vx = func_eval(val_x,x_i,y_j);
        vy = func_eval(val_y,x_i,y_j);
        %setting up the vector with actual solutions!!
        %solving the given equation subject to the boundary conditions that
        %u(x,y) = 0 at all the boundaries, we get u(x,y) =
        %sin(pi*x)sin(pi*y)/(2pi^2)
        
        u_temp(i,j) = (sin(pi*x_i)*sin(pi*y_j))/(2*(pi^2));
        u_actual(count)= u_temp(i,j);          % converting the u matrix into a vector
        f_actual(count)= func_eval(f,x_i,y_j); % to generate the RHS of Au = f
        
        % setting up co-efficients for the nodes at the edges
         if(i==1 && j==1) %point (1,1)
            A(count,count) = 1;
            
        elseif (i==1 && j==m) %point(1,m)
           A(count,count) = 1;
        elseif (i==m && j==1) %point(m,1)
            A(count,count) = 1;
            
         elseif (i==m && j==m) %point(m,m)
            A(count,count) = 1;
            
        % end of setting up co-efficients for nodes
        
        %setting up coefficients for points ON the boundary line
               
        elseif (i==1 && i~=j && j~=m)
           A(count,count) = 1;
            
        elseif (i==m && i~=j && j~=1)
            A(count,count) = 1;
        elseif(j==1 && i~=j && i~=m)
            A(count,count) = 1;
            
        elseif(j==m && i~=j && i~=1)
           A(count,count) = 1;
           
        else %setting up the co-efficients for the points in the interior
            A(count,count) = 4/(h^2);
            A(count,(count-1)) = -1/(h^2) - (vx/(2*h));
            A(count,(count+1)) = -1/(h^2) + (vx/(2*h));
            A(count,(count-m)) = -1/(h^2) - (vy/(2*h));
            A(count,(count+m)) = -1/(h^2) + (vy/(2*h));
                   
        end
     end
end
sol1 = A\f_actual;           %solving the A=f system using the intrinsic solver
calc_grid = reshape(sol1,m,m);% converting the solution vector to a matrix so that a contour plot can be generated
sol2 = calc_grid';           
contour(linspace(0,1,m),linspace(0,1,m),sol2);
calc_err = u_actual - sol1;
sol3 = norm(calc_err,2); %Calculating the L2 norm of error
end

