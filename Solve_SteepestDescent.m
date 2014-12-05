function [ X ] = Solve_SteepestDescent ( A, b, x0, e )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    X = x0';
    r = b - A*X;
    while norm(r, inf)>e
        X = X + (r'*r)/(r'*A*r)*r;
        r = b - A*X
    end
end

