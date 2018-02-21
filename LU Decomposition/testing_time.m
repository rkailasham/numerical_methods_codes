function testing_time

tic
disp('First call');
A=rand(5);
B=rand(5);
A*B;
a=toc;
tic
disp('second call');
b = toc;
disp(a);
disp(b);
end