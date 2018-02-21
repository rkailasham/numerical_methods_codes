%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A - The LU decomposed version of the A in Ax = b                         %
%O - Array of size [1 x n] that contains info about row interchanges      %
%n - dimension of A                                                       %
%b - RHS of the equation Ax = b                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% OUTPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sol - A [1 x n] array that contains the solution x of Ax=b               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol] = Substitute (A,O,n,b)
x = zeros(1,n);
%This part of the code takes care of forward substitution
for i =2:n
    sum = b(O(i));
    for j = 1 : (i-1)
        sum = sum - A(O(i),j)*b(O(j));
    end
    b(O(i)) = sum;
end
%End of forward substitution routine
x(n) = b(O(n))/A(O(n),n);

%This part does the backward substitution
for i = (n-1) : -1 : 1
    sum = 0;
    for j = (i+1) : n
        sum = sum + A(O(i),j)*x(j);
    end
    x(i) = (b(O(i)) - sum)/A(O(i),i);
end
%End of backward substitution
sol = x;

end


    
