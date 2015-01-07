function [ E, I ] = Solve_Eig_Givens( A, e )
% 使用上 Hessenberg 分解, Givens 变化求解 A 的特征值
%[ E, I ] = Solve_Eig_Givens( A )
%   A 方阵
%   e 精度
%返回值
%   E 特征值列向量
%   I 特征值对角阵(与E列对应)

n = length(A);
T = A;
for k = 1:n-2
    [v,b] = Solve_Householder(T(k+1:n,k));
    H = eye(n-k) - b*v'*v;
    T(k+1:n,k:n) = H*T(k+1:n,k:n);
    T(1:n,k+1:n) = T(1:n,k+1:n)*H;
end
C = Get_Tril(T,n);
T = Special_Givens(T,n);
c = Get_Tril(T,n);
k = 0;

while norm(C-c,Inf)>e
    k = k + 1;
    if k>5000; break; end
    C = c;
    T = Special_Givens(T,n);
    c = Get_Tril(T,n);
end
k = 1;
E = zeros(n);
I = zeros(n);
while k<n
    if abs(T(k+1,k)) > e
        b = T(k,k)+T(k+1,k+1);
        I(k,k) = (b+sqrt(b^2-4*(T(k,k)*T(k+1,k+1)-T(k,k+1)*T(k+1,k))))/2;
        I(k+1,k+1) = conj(I(k,k));
        [l,E(:,k)] = Solve_BackwardPowerMethod(A,I(k,k)-e,e);
        E(:,k+1) = conj(E(:,k));
        k = k + 2;
    else
        [I(k,k),E(:,k)] = Solve_BackwardPowerMethod(A,T(k,k)-e,e);
        k = k + 1;
    end
end
if k==n; [I(k,k),E(:,k)] = Solve_BackwardPowerMethod(A,T(k,k)-e,e); end
end

function [ B ] = Special_Givens( A, n )
    B = eye(n);
    for i = 1:n-1
        [c,s] = Solve_Givens(A(i,i), A(i+1,i));
        D = eye(n);
        D(i,i) = c;
        D(i,i+1) = s;
        D(i+1,i) = -s;
        D(i+1,i+1) = c;
        B = D*B;
        A = D*A;
    end
    B = A*B';
end

function [c] = Get_Tril( A, n )
    c = zeros(n-1,1);
    for i = 1:n-1
        c(i) = A(i+1,i);
    end
end

function tester()
    n = 5;
    A = rand(n);
    [E,I] = Solve_Eig_Givens(A,0.0000001);
    diag(I)
    for i = 1:n
        disp(norm(A*E(:,i)-I(i,i)*E(:,i)));
    end
end
