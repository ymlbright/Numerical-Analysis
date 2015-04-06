function [X] = Solve_Chasing(d, e, f, b)
%使用追赶法求解方程组 Ax=b
% d     系数矩阵次对角线系数(列向量)
% e     系数矩阵对角线系数(列向量)
% f     系数矩阵超对角线系数(列向量)
% b     右端向量
%返回值
% X     方程组的解 
n = length(b);
X = zeros(n,1);
for i = 2:n
    l = d(i-1)/e(i-1);
    e(i) = e(i) - f(i-1)*l;
    b(i) = b(i) - b(i-1)*l;
end
X(n) = b(n)/e(n);
for i=n-1:-1:1
    X(i) = (b(i)-f(i)*X(i+1))/e(i);
end
end