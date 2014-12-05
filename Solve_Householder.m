function [v, b] = Solve_Householder(x)
%寻找 Householder 矩阵 H 使 Hx 的后 n-1 个分量为 0
%[v, b] = Solve_Householder(x)
%   x  列向量
%返回值:
%   [v,b] H=I-b*v'*v
    n = length(x);
    x = x/norm(x, inf);
    s = x(2:n)'*x(2:n);
    v(2:n) = x(2:n);
    if s == 0
        b = 0;
    else
        a = sqrt(x(1)^2 + s);
        if x(1) <= 0
            v(1) = x(1) - a;
        else
            v(1) = -s/(x(1)+a);
        end
        b = 2*v(1)^2/(s+v(1)^2);
        v = v/v(1);
    end
end

