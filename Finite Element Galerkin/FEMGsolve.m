%***********************************************************************%
% This function solves the given one-dimensional differential equation  %
% using Finite Element Method with linear Galerkin elements.            %
% The function takes the number of elements, the Diffusivity, the       %
% velocity, the rate constant and the tolerance as input. It returns    %
% c(x) as the solution.                                                 %
%                                                                       %
%***********************************************************************%


function [sol] = FEMGsolve(N_el,D,v,k,tol)

N = N_el+1;
h =1/N_el;
syms e;

x = zeros(1,N);
for i=1:N
    x(i) = (i-1)*h;
end

c = 0.5*ones(N,1);
%Mapping between global and master element
phi(1) = 0.5*(1-e);
phi(2) = 0.5*(1+e);

dphi_de(1) = diff(phi(1),e);
dphi_de(2) = diff(phi(2),e);

norm_delc = 1;
count = 0;

c(1) = 1;  % boundary condition. Has to be true.
while(norm_delc>tol) % Beginning of Newton's loop
    J = zeros(N,N);
    R = zeros(N,1);
    J_loc = zeros(2,2);
    R_loc = zeros(2,1);
    x_loc = zeros(2,1);
    
    for l=1:N_el % looping over all elements
      T = sym(zeros(2,1));  % for evaluating residual
      S = sym(zeros(2,2)); % for evaluating Jacobian      
      c_loc(1) = c(l);
      c_loc(2) = c(l+1);
      x_loc(1) = x(l);
      x_loc(2) = x(l+1);
      
      Nlin =sym(zeros(2,1));
      
      X = x_loc(1)*dphi_de(1) + x_loc(2)*dphi_de(2);
      for i=1:2
          
            for j=1:2
                
            
            T(i) = T(i) + c_loc(j)*(D*(dphi_de(i)*dphi_de(j)/X) +v*(phi(i)*dphi_de(j)));
            if(j==1)
                flag=2;
            else flag=1;
            end
            S(i,j) = D*(dphi_de(i)*dphi_de(j)/X) +v*(phi(i)*dphi_de(j)) + k*phi(i)*X* (2*c_loc(j)*phi(j)^2 + 2*c_loc(flag)*phi(1)*phi(2));
            J_loc(i,j) = GaussEval(S(i,j));
            end
            
            Nlin(i) = phi(i)*k* X*(c_loc(1)^2*phi(1)^2 + 2*c_loc(1)*c_loc(2)*phi(1)*phi(2) + c_loc(2)^2*phi(2)^2);
            R_loc(i) = GaussEval(T(i) + Nlin(i));
        end
        R(l) = R(l) + R_loc(1);
        R(l+1) = R_loc(2);
        
        J(l,l) = J(l,l) + J_loc(1,1);
        J(l,l+1) = J_loc(1,2);
        J(l+1,l) = J_loc(2,1);
        J(l+1,l+1) = J_loc(2,2);
    end
        
   R(1) = 0;
   J(1,1) = 1;
   J(1,2) = 0;
   R = -1*R;
   delc = J\R;
   c = c+delc;
   norm_delc = norm(delc);
   count=count+1;         
   disp(norm_delc);   
   end
sol = c;

disp('Done with one solve');
end

