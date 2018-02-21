%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to test the performance of user-written LU code on the system  %
% (A^n)x = b. A an b are randomly generated at the beginning. The exponent%
% 'n' is continually varied and the behavior of x is observed.            %
%  A system of 25 equations is studied.                                   %
% The function returns the following :                                    %
% x    - A [25 x n] array that contains the solution to (A^n)x = b for    %
%        each value of n.
% time - Time taken by the user-written code to solve the system(A^n)x = b%
%        for each value of n.                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [x, time] = Test_HigherPower(n)
x = zeros(25,n);
time = zeros(1,n);
A = randi(12,25);
b = [1:25]';
B = eye(25);       %creates an identity matrix of size 25 x 25
for i=1:n
    B=B*A;         %effectively, B = A^n
    tic            %begin stopwatch to track Kailash's LU code 
    x(:,i)= Ludecomp(B,b,25,10e-30);
    time(i) = toc; %end stopwatch to track Kailash's LU code
end
end

    