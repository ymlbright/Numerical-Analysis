function [L, U, P] = Copy_of_Decompose_LU(A)
%使用Gauss法进行LU矩阵分解
%   A 方阵
%[L, U, P] = Decompose_LU(A)
%返回值:
%   L 下三角矩阵
%   U 上三角矩阵,其对角线上为 D
%   P 置换矩阵
n = length(A);
P = eye(n);
L = zeros(n);
U = zeros(n);
for k = 1:n-1
    if A(k,k) ~= 0
        A(k+1:n,k) = A(k+1:n,k)/A(k,k);
        A(k+1:n,k+1:n) = A(k+1:n,k+1:n) - A(k+1:n,k)*A(k,k+1:n);
    else
        disp('矩阵A奇异！')
        pause
    end
end
for i = 1:n
    for j = 1:n
        if i>j
            L(i,j) = A(i,j);
        else
            U(i,j) = A(i,j);
        end
    end
end
L = L + eye(n);
end
