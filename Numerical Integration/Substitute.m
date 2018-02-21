%A - Co-efficient array
%O - Array of size 1 x n that keeps track of row interchanges
%b - RHS of the equation Ax = b
%x - vector containing unknowns x1,x2,...xn
% x and b have the same dimensions.
%n - dimension of A
%k - passed in as a parameter. Part of the bigger picture.

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


    
