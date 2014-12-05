function [ X ] = Solve_QR( A, b )
%使用QR分解法求解超定线性方程组
%[ X ] = Solve_QR( A, b )
%   A  超定方程组系数矩阵
%   b  方程组右端值(列向量)
%返回值:
%   X 方程组的解
[Q, R] = Decompose_QR(A);
X = Solve_U(R, Q'*b);

end

