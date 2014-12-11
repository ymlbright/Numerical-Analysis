function [X] = Copy_of_Solve_Gauss(A, b)
%使用高斯消去法求解方程组 Ax=b
%[X] = Copy_of_Solve_Gauss(A, b)
%   A 方程组系数矩阵
%   b 方程组右端值
%返回值:
%   X 方程组的解
[L, U, P] = Copy_of_Decompose_LU(A);
X = Solve_L(L, P*b);
X = Solve_U(U, X);
end