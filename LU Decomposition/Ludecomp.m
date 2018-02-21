%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A - Co-efficient array i.e the A in Ax = b                               %
%b - RHS of the equation Ax = b                                           %
%n - dimension of A                                                       %
%tol - Tolerance value as decided by the user. 10e-6 will suffice for     %  
%      normal operations and commonly encountered matrices (lower tol     %
%      value(10e-30) needed for herbert matrices.                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% OUTPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sol - A [1 x n] array that contains the solution x of Ax=b               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [sol] = Ludecomp(A,b,n,tol)
O = zeros(1,n);   %Array of size 1 x n that keeps track of row interchanges
S = zeros(1,n);   %Array of size 1 x n that keeps track of the largest element in each row
er = 0;
[er,O,S,A] = Decompose(A,n,tol,O,S);
if er~=-1
    sol = Substitute(A,O,n,b);
else disp('Tolerance limit exceeded');
     sol = zeros(1,n);
end
end
