function [ c, s ] = Solve_Givens( a, b )
%Givens±ä»»
%[ c, s ] = Solve_Givens( a, b )
    if b == 0
        c = 1; s = 0;
    else
        if norm(b, inf) > norm(a, inf)
            f = a/b;
            s = 1/sqrt(1 + f^2);
            c = s*f;
        else
            f = b/a;
            c = 1/sqrt(1 + f^2);
            s = c*f;
        end
    end
end

