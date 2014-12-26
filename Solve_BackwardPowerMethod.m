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

Pn = 10;
Pool = zeros(n,Pn);
Pool(:,Pn) = Solve_L(L, P*x);
Pool(:,Pn) = Solve_U(U, Pool(:,Pn));
exit_flag = 0;

while norm(Pool(:,Pn)-Pool(:,1))>e && ~exit_flag
    if k > 10000; disp('超过最大迭代次数,反幂法可能不收敛.'); break; end
    Pool(:,1) = Pool(:,Pn);
    for i = 2:Pn
        Pool(:,i) = find_next(Pool(:,i-1), L, U, P, n);
        if norm(Pool(:,i)-x)<e; exit_flag = 1;break; end
    end
    k = k + Pn;
end

X = Pool(:,Pn);
m = -Inf;
for i = 1:n ; if m<X(i); m = X(i);j=i; end; end;
lamda = A(j,:)*X/X(j);
X = X./norm(X);

end

function [X]=find_next(x, L, U, P, n)
    y = Solve_L(L, P*x);
    y = Solve_U(U, y);
    m = -Inf;
    for i = 1:n ; if m<y(i); m = y(i); end; end;
    X = y/m;
end