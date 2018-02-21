%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A - Co-efficient array i.e the A in Ax = b                               %
%O - Array of size [1 x n] that keeps track of row interchanges           %
%S - Array of size [1 x n] that keeps track of the largest element in each% 
%    row                                                                  %
%n - dimension of A                                                       %
%k - passed in as a parameter. Keeps track of the row which is being      %
%     checked for pivoting.                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% OUTPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sol - A [1 x n] array that contains information about the order of rows %
%        after pivoting.                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol] = Pivot(A,O,S,n,k)

p=k;
big = abs(A(O(k),k)/S(O(k))); % each row element scaled by didviding by the 
                              %  largest element in the row.
for ii = (k+1) : n
    dummy = abs(A(O(ii),k)/S(O(ii)));
    if dummy > big   %if the scaled value in a subsequent row is bigger than
        big = dummy; % the scaled value of the previous row, then pivot! 
        p = ii;
    end 
end 
dummy = O(p);
O(p) = O(k);
O(k) = dummy;
sol = O;
end 