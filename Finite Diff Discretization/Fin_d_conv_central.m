function [sol1, sol2, sol3] = Fin_d_conv(n,val_x,val_y,f,ref_true,mode)
h = 1/n;
N = (n+1)^2;
m = (n+1);
A = zeros(N,N);
factor = 100/m;
x_i=0;
y_j=0;
syms x;
syms y;
vx=0;
vy=0;
%A(1,1) = 1;
%A((n+1),(n+1)) = 1;
l = zeros (1,N);
u = zeros (N,1);
%u_temp = zeros(n+1);
%u_actual = zeros(N,1);
f_actual = zeros(N,1);
calc_error = zeros(m,m);
count=0;
i_ext = 0;
j_ext = 0;
for i=1:m
    for j=1:m
        %i_ext = ((i-1)*(101-1)/(N-1)) + 1;
        %j_ext = ((i-1)*(101-1)/(N-1)) + 1;
        count = ((j-1)*(n+1) + i);
        x_i = (i-1)*h;  %setting up the vector with actual solutions!!
                        %solving the given equation subject to the boundary conditions that
                        %u(x,y) = 0 at all the boundaries, we get u(x,y) =
                        %sin(pi*x)sin(pi*y)/(2pi^2)
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
        %using central difference for the convection term (first
        %derivative) and appropriate method for second derivative.
        
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
sol1 = A\f_actual;
calc_grid = reshape(sol1,m,m);
sol2 = calc_grid';
for i=1:m
    for j=1:m
      i_ext = (((i-1)*(151-1)/(m-1))) + 1;
      j_ext = (((j-1)*(151-1)/(m-1))) + 1; 
      calc_error(i,j) = sol2(i,j) - ref_true(i_ext,j_ext);
    end
end

err_vector = reshape(calc_error,N,1);
sol3 = (norm(err_vector,2)/m);
end

