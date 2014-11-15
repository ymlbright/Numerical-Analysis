function [r] = Estimate_Reverse_Matrix_Modinf(A, m)
%估计矩阵 A 的逆的无穷范数
%[r] = Estimate_Reverse_Matrix_Modinf(A, m)
% A     n*n方阵
% m     精度,最大逼近次数
%返回值
% r     矩阵 A 的逆的无穷范数估计值
r = Estimate_Reverse_Matrix_Mod1(A', m);
end