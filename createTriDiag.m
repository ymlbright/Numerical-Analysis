function [A] = createTriDiag(d, e, f, n)
%构造三对角矩阵A
%[A] = createTriDiag(d, e, f, n)
% d     系数矩阵次对角线系数(列向量)
% e     系数矩阵对角线系数(列向量)
% f     系数矩阵超对角线系数(列向量)
% n     A的大小
%返回值
% A     三对角矩阵
    A = zeros(n);
    for i = 1:n
        A(i,i) = e(i);
        if i<n ; A(i,i+1) = f(i); end
        if i>1 ; A(i,i-1) = d(i-1); end
    end
end

