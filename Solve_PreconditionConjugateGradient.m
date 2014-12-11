function [ X ] = Solve_PreconditionConjugateGradient( A, b, M, x0, e )
%使用预优共轭梯度法求解线性方程组
%[ X ] = Solve_PreconditionConjugateGradient( A, b, M, x0, e )
%   A  对称正定矩阵
%   b  方程组右端值(列向量)
%   M  预优矩阵
%   x0 初始迭代点(行向量)
%   e  误差限
%返回值:
%   X 方程组的解
X = x0';
r = b - A*X;
k = 0;
while sqrt(r'*r)>e
    z = Solve_Gauss(M, r);
    k = k + 1;
    if k > 100000; disp('ConjugateGradient 超过最大迭代次数!'); break; end
    if k == 1
        p = z;
        rou1 = r'*z;
    else
        rou2 = rou1;
        rou1 = r'*z;
        beta = rou1/rou2;
        p = z + beta*p;
    end
    w = A*p;
    alpha = rou1/(p'*w);
    X = X + alpha*p;
    r = r - alpha*w;
end
end

