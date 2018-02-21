
%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A - Co-efficient array i.e the A in Ax = b                               %
%O - Array of size [1 x n] that keeps track of row interchanges           %
%S - Array of size [1 x n] that keeps track of the largest element in each% 
%    row                                                                  %
%n - dimension of A                                                       %
%tol - tolerance value. As entered by the user.                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% OUTPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%flag - set equal to 0 if tolerance value is not violated. -1 otherwise   %
%Order- Contains the same information as in O.                            %  
%Si   - Contains the same information as in S.                            %
%A_Changed - Matrix A after it has been LU Decomposed                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [flag, Order, Si, A_changed] = Decompose(A,n,tol,O,S)
flag=0;
for i=1:n
    O(i) = i;
    S(i) = abs(A(i,1));
    for j=2:n
        if abs(A(i,j)) > S(i)
            S(i) = abs(A(i,j)); %stores the largest element in each row
        end
    end
end
for k=1:(n-1)
    O = Pivot(A,O,S,n,k); % Keeps track of row-interchanges during pivoting
    if abs(A(O(k),k)/S(O(k))) < tol
       flag = -1;
       disp(A(O(k),k)/S(O(k)));
       Order = O;
       Si = S;
       A_changed = A;
       return;
    end
    for i= (k+1) : n
        factor = A(O(i),k)/A(O(k),k);
        A(O(i),k) = factor;
        for j = (k+1) : n
            A(O(i),j)=A(O(i),j) - factor*A(O(k),j);
        end
    end
end
if abs(A(O(k),k)/S(O(k))) < tol
    flag = -1;
    disp(A(O(k),k)/S(O(k)));
else flag = 0;
end
Order = O;
Si = S;
A_changed = A;
end


