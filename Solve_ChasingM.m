function [X] = Solve_ChasingM(D, E, F, m, b)
%使用追赶法求解方程组
%| E F   |
%| D E F | x = b
%|   D E |
% D     次对角线矩阵
% E     主对角线矩阵
% F     超对角线矩阵
% m     主对角线矩阵重复次数
% b     右端向量
%返回值
% X     方程组的解

% [
%     2 3 0 1 0 0 0 0 0;
%     3 2 3 0 1 0 0 0 0;
%     0 3 2 0 0 1 0 0 0;
%     1 0 0 2 3 0 1 0 0;
%     0 1 0 3 2 3 0 1 0;
%     0 0 1 0 3 2 0 0 1;
%     0 0 0 1 0 0 2 3 0;
%     0 0 0 0 1 0 3 2 3;
%     0 0 0 0 0 1 0 3 2;
% ]

    n = length(E);
    if m~=n; disp('只支持方阵计算！'); X=0; return; end
    if m == 1
        X = [E F;D E]\b;
    else
        A = zeros(m*n,m*n);
        A(1:n,1:2*n) = [E F];
        for i=1:m-2; A(1+i*n:(i+1)*n,1+(i-1)*n:(i+2)*n) = [D E F]; end
        A(1+(m-1)*n:m*n,1+(m-2)*n:m*n) = [D E];
        X = A\b;
    end
end