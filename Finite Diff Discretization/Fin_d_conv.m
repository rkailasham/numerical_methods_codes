%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to solve the system Au = f for the case with      %
% convection using the central or forward difference schemes for the      %
% convection term as chosen by the user.                                  %
%                                                                         %
%%%%%%%%%%%%%%%% INPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% n        : Number of grid points to be used                             %
% val_x    : The x-component of velocity.                                 %
% val_y    : The y-component of velocity.                                 %
% f        : The forcing function                                         %
% ref_true : The reference solution with which the results are to be      %
%            compared                                                     %
% ref_N    : The size of the reference grid. It is 151 for the central    %
%            difference case and 121 for the forward difference case      %
% mode     : Can be either 'central' or 'forward'.                        %
%                                                                         %
%%%%%%%%%%%%%%%%% FUNCTION OUTPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% sol1     : Contains the solution u of the system Au = f                 % 
% sol2     : Contains the matrix version of the solution u.               %
% sol3     : Contains the solution error normalised the number of grid    %
%            points (i.e (n+1) )                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol1, sol2, sol3] = Fin_d_conv(n,val_x,val_y,f,ref_true,ref_N,mode)
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
f_actual = zeros(N,1);
calc_error = zeros(m,m);
count=0;
i_ext = 0;
j_ext = 0;
for i=1:m
    for j=1:m
        count = ((j-1)*(n+1) + i); %Converting the double-index system to a single index
        x_i = (i-1)*h;  
        y_j = (j-1)*h;
        vx = func_eval(val_x,x_i,y_j);
        vy = func_eval(val_y,x_i,y_j);
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
           
        else 
            if(strcmp(mode,'central')
                A(count,count) = 4/(h^2);
                A(count,(count-1)) = -1/(h^2) - (vx/(2*h));
                A(count,(count+1)) = -1/(h^2) + (vx/(2*h));
                A(count,(count-m)) = -1/(h^2) - (vy/(2*h));
                A(count,(count+m)) = -1/(h^2) + (vy/(2*h));
            
            %forward difference
            elseif(strcmp(mode,'forward'))
                A(count,count) = 4/(h^2) - vx/h - vy/h;
                A(count,(count-1)) = -1/(h^2);
                A(count,(count+1)) = -(1/(h^2)) + (vx/h);
                A(count,(count-m)) = -1/(h^2);
                A(count,(count+m)) = -(1/(h^2)) + (vy/h);
            end
        end
     end
end
sol1 = A\f_actual;             %solving the A=f system using the intrinsic solver
calc_grid = reshape(sol1,m,m); % converting the solution vector into a matrix for plotting
sol2 = calc_grid';       
for i=1:m
    for j=1:m
      i_ext = (((i-1)*(ref_N-1)/(m-1))) + 1; %linear extrapolation to map
      j_ext = (((j-1)*(ref_N-1)/(m-1))) + 1; %points from the coarse matrix to the fine matrix
      calc_error(i,j) = sol2(i,j) - ref_true(i_ext,j_ext);
    end
end

err_vector = reshape(calc_error,N,1);
sol3 = (norm(err_vector,2)/m);  %Calculating the L2 norm of error
end

