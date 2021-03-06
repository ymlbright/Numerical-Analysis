function [X] = Solve_U(U, b)
%使用回代法求解上三角形方程组
%[X] = Solve_U(U, b)
%   U 上三角形矩阵
%   b 方程组右端值
%返回值:
%   X 方程组的解
[m,n] = size(U);
n = min([m n]);
for j = n:-1:2
    if U(j,j)==0
        b(j) = 0;
    else
        b(j) = b(j)/U(j,j);
    end
    b(1:j-1) = b(1:j-1) - b(j)*U(1:j-1,j);
end
b(1) = b(1)/U(1,1);
X = b(1:n);
end