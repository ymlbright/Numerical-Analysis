function [ X ] = Solve_SOR( A, b, x0, w, e )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    n = length(A);
    [L, D, U] = ToLDU(A, n);
    M = (D-w*L)\((1-w)*D+w*U);
    g = w*((D-w*L)\b);
    x = x0';
    X = M*x + g;
    while norm(X-x, inf)>e
        x = X;
        X = M*x + g
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
