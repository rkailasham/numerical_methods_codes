function [ret] = time_track
t=zeros(1,100);

for i=1:100
    tic
    A=randi((i*10),i);
    B=randi((i*5),i);
    A*B;
    t(i)=toc;
end
ret=t;
end
