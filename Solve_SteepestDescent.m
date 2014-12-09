function [ X ] = Solve_SteepestDescent ( A, b, x0, e )
%使用最速下降法求解线性方程组
%[ X ] = Solve_SteepestDescent ( A, b, x0, e )
%   A  对称正定矩阵
%   b  方程组右端值(列向量)
%   x0 初始迭代点(行向量)
%   e  误差限
%返回值:
%   X 方程组的解
    X = x0';
    r = b - A*X;
    while norm(r, inf)>e
        X = X + (r'*r)/(r'*A*r)*r;
        r = b - A*X;
    end
end

