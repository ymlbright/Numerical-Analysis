function [X] = Solve_L(L, b)
%使用前代法求解下三角形方程组
%[X] = Solve_L(L, b)
%   L 下三角形矩阵
%   b 方程组右端值
%返回值:
%   X 方程组的解
n = length(L);
for j = 1:n-1
    b(j) = b(j)/L(j,j);
    b(j+1:n) = b(j+1:n) - b(j)*L(j+1:n,j);
end
b(n) = b(n)/L(n,n);
X = b;
end