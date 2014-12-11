function [ lamda, X ] = Solve_BackwardPowerMethod( A, u, e )
%使用反幂法求 A 中特征值距离 u 最近的特征值及对应向量
%   A 需要求解的矩阵
%   u 给定值
%   e 精度
%返回值
%   lamda   特征值
%   X       特征值对应的向量

n = length(A);
x = rand(n,1);
x(1) = 1; j = 1; k = 0;

% 进行LU分解
[L, U, P] = Decompose_LU(A-u*eye(n));

X = Solve_L(L, P*x);
X = Solve_U(U, X);

while norm(X-x)>e
    if k > 10000; disp('超过最大迭代次数,反幂法不可能收敛.'); break; end
    k = k + 1;
    y = Solve_L(L, P*X);
    y = Solve_U(U, y);
    m = -Inf;
    for i = 1:n ; if m<y(i); m = y(i);j = i; end; end;
    x = X;
    X = y/m;
end
lamda = A(j,:)*X/X(j);
end

