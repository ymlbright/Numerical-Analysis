function [X] = Copy_of_Solve_Cholesky(A, b)
%使用平方根法求解方程组 Ax=b
%[X] = Solve_Cholesky(A, b)
%   A 方程组系数矩阵
%   b 方程组右端值
%返回值:
%   X 方程组解的解
L = Copy_of_Decompose_Cholesky(A);
X = Solve_L(L, b);
X = Solve_U(L', X);
end