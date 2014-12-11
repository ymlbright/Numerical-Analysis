function [A] = createHilbert(n)
%¹¹Ôìn½×Hilbert¾ØÕó
%[A] = createHilbert(n)
    A = zeros(n,n);
    for i = 1:n
        for j  = i:n
            A(i,j) = 1/(i+j-1);
            A(j,i) = A(i,j);
        end
    end
end