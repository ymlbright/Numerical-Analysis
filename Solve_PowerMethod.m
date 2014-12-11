function [ lamda, X ] = Solve_PowerMethod( A, e )
%使用幂法求 A 最大特征值
%   A 需要求解的矩阵
%   e 精度
%返回值
%   lamda   最大特征值
%   X       最大特征值对应的向量

n = length(A);
x = rand(n,1);
x(1) = 1; j = 1; k = 0;
X = A*x;
while norm(X-x)>e
    if k > 10000; disp('超过最大迭代次数,幂法不可能收敛.'); break; end
    k = k + 1;
    y = A*X;
    m = -Inf;
    for i = 1:n ; if m<y(i); m = y(i);j = i; end; end;
    x = X;
    X = y/m;
end
lamda = A(j,:)*X/X(j);
end

