function [ c, s ] = Solve_Givens( a, b )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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

