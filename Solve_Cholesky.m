function [X] = Solve_Cholesky(A, b)
%使用改进的平方根法求解方程组 Ax=b
%[X] = Solve_Cholesky(A, b)
%   A 方程组系数矩阵
%   b 方程组右端值
%返回值:
%   X 方程组解的解
[L, D] = Decompose_Cholesky(A);
X = Solve_L(L, b);
X = Solve_U(D*L', X);
end