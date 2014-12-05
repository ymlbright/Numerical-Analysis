function [ X ] = Solve_QR( A, b )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[Q, R] = Decompose_QR(A);
X = Solve_U(R, Q'*b);

end

