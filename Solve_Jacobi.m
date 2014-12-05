function [ X ] = Solve_Jacobi( A, b, x0, e )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    n = length(A);
    [L, D, U] = ToLDU(A, n);
    M = D\(L+U);
    g = D\b;
    x = x0';
    X = M*x + g;
    while norm(X-x, inf)>e
        x = X;
        X = M*x + g;
    end
end

function [L, D, U] = ToLDU(A, n)
    L = zeros(n);
    D = zeros(n);
    U = zeros(n);
    for i = 1:n
        L(i,1:i-1) = -A(i,1:i-1);
        D(i,i) = A(i,i);
        U(i,i+1:n) = -A(i,i+1:n);
    end
end