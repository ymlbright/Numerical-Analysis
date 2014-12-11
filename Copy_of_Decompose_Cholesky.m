function [L] = Copy_of_Decompose_Cholesky(A)
%使用平方根法对矩阵A进行分解
%[L, U] = Decompose_Cholesky(A)
%   A 方阵
%返回值:
%   L 下三角阵
%   D 对角线元素                 
n = length(A);
for k = 1:n
    A(k,k) = sqrt(A(k,k));
    A(k+1:n,k) = A(k+1:n,k)/A(k,k);
    for j = k+1:n
        A(j:n,j) = A(j:n,j) - A(j:n,k)*A(j,k);
    end
end
L = tril(A);
end