function [ E, I ] = Solve_Eig( A, e )
% 求解 A 的特征值
%[ E, I ] = Solve_Eig( A )
%   A 方阵
%   e 精度
%返回值
%   E 特征值列向量
%   I 特征值对角阵(与E列对应)
if A' == A
    [E,I] = Solve_Eig_Givens_Symmetric(A,e);
else
    [E,I] = Solve_Eig_Givens(A,e);
end
end