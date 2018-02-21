%Function to generate Hilbert matrix of size n x n.
%Elements of a Hilbert matrix obey the rule H(i,j) = 1/(i+j-1)

function [H] = hilbert(n)
h=zeros(n,n);
for i =1:n
    for j=1:n
       h(i,j) = 1/(i+j-1);
    end
end
H = h;
end
