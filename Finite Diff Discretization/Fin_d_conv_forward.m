function [sol1, sol2,sol3,sol4] = Fin_d_conv_forward(n,val_x,val_y,f,ref_true)
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
calc_error = zeros(N,N);
count=0;
for i=1:m
    for j=1:m
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
            A(count,count) = 4/(h^2) - vx/h - vy/h;
            A(count,(count-1)) = -1/(h^2);
            A(count,(count+1)) = -(1/(h^2)) + (vx/h);
            A(count,(count-m)) = -1/(h^2);
            A(count,(count+m)) = -(1/(h^2)) + (vy/h);
                   
        end
     end
end
sol1 = A\f_actual;
calc_grid = reshape(sol1,m,m);
sol2 = calc_grid';
%contourf(linspace(0,1,m),linspace(0,1,m),sol2);

x_i=0;
y_j=0;
X_k=0;
Y_p=0;
for i=1:m
    for j=1:m
        x_i = (i-1)*h;
        y_j = (j-1)*h;
        for k=1:101
            for p=1:101
                X_k = (k-1) * 0.01;
                Y_p = (p-1) * 0.01;
        if((x_i == X_k) && (y_j == Y_p))
            calc_grid(i,j) = ref_true(k,p) - sol2(i,j);
        end
            end
        end
    end
end
re_grid = reshape(calc_grid,N,1);
sol3 = norm(re_grid,2);
sol4 = calc_grid;

end

