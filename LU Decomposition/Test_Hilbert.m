%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to test the performance of the inbuilt and user-defined        %
% LU codes on Hilbert matrices of different sizes.                        %
% Takes the maximum size of matrix that needs to be tested as an argument %
% and returns the following :                                             %
%
% C_No      - Condition number of the matrix being solved                 %
% time_user - Time taken by the user-written code to solve the system Hx=b%
% time_comp - Time taken by the inbuilt LU function to solve the same 
%             system Hx = b                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [C_No, time_user, time_comp] = Test_Hilbert(n)
time_comp = zeros(1,n);
time_user = zeros(1,n);
S = zeros(1,n);
O = zeros(1,n);
for i=1:n
    H = hilbert(i);
    b= [1:i]';
    [L,U,P]=lu(H);
    o = [1:i];
    tic                 %begin stopwatch to track computer's LU code
    B = lu(H);
    time_comp(i) = toc; %end stopwatch to track computer's LU code
    
    tic                 %begin stopwatch to track Kailash's LU code
    [flag,O,S,A] = Decompose(H,i,10e-50,O,S);
    time_user(i) = toc; %end stopwatch to track Kailash's LU code
    C_No(i)=cond(H);
    
end
end
