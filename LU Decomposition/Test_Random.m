%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to test the performance of the inbuilt and user-defined        %
% LU codes on Random matrices of different sizes.                         %
% Takes the maximum size of matrix that needs to be tested as an argument %
% and returns the following :                                             %
%
% C_No      - Condition number of the matrix being solved                 %
% time_user - Time taken by the user-written code to solve the system Rx=b%
% time_comp - Time taken by the inbuilt LU function to solve the same 
%             system Rx = b                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function [C_No, time_user, time_comp] = Test_Random(n)
time_comp = zeros(1,n);
time_user = zeros(1,n);
S = zeros(1,n);
O = zeros(1,n);
flag = 0;
for i=1:n
    flag = 0;
    R = randi((4*i),i);
    b= [1:i]';
    o = [1:i];
    %tic
    tic                 %begin stopwatch to track computer's LU code
    B = lu(R);
    time_comp(i) = toc; %end stopwatch to track computer's LU code
    
    tic                 %begin stopwatch to track Kailash's LU code 
    [flag,O,S,A] = Decompose(R,i,10e-50,O,S);
    time_user(i) = toc; %end stopwatch to track Kailash's LU code
    C_No(i)=cond(R);
    
end
end
 
