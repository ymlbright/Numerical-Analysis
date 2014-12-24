function task3()
    clc
    %fprintf(1,'N=20\t\t%f\n',norm(p137_2(20)-1));
    %fprintf(1,'N=40\t\t%f\n',norm(p137_2(40)-1));
    %fprintf(1,'N=80\t\t%f\n',norm(p137_2(80)-1));
    
%     fprintf(1,'N=20\t\t\n');
%     p159_1(20);
%     fprintf(1,'N=40\t\t\n');
%     p159_1(40);
%     fprintf(1,'N=80\t\t\n');
%     p159_1(80);

    
%     fprintf(1,'N=20\t\t\n');
%     p159_1_0(20);
%     fprintf(1,'N=40\t\t\n');
%     p159_1_0(40);
%     fprintf(1,'N=80\t\t\n');
%     p159_1_0(80);

%     fprintf(1,'N=20\t\t%f\n',p160_2(20));
%     fprintf(1,'N=40\t\t%f\n',p160_2(40));
%     fprintf(1,'N=80\t\t%f\n',p160_2(80));

    p160_3
end

function [r]=f137_2(x,y)
    r=x+y;
end

function [r]=g137_2(x,y)
    r=exp(x*y);
end

function [A] = createA(S,B,n)
% 生成n阶三对角对称矩阵, 次对角线为矩阵B, 主对角线矩阵由函数 S 生成
    m = length(B);
    A = zeros(m*n);
    for i = 0:n-1
        if i>0 ; A((i-1)*m+1:i*m,i*m+1:(i+1)*m)=B; end
        A(i*m+1:(i+1)*m,i*m+1:(i+1)*m)=S(i,m); 
        if i<n-1 ; A((i+1)*m+1:(i+2)*m,i*m+1:(i+1)*m)=B; end
    end
end

function [b] = b137_2(i,j)
% 计算 P137(2) 差分方程组右端向量
    global h_137_2; 
    b = h_137_2^2 * f137_2(i*h_137_2, j*h_137_2);
end

function [u] = fu137_2(i,j)
    u = 1;
end

function [b] = createb(fb, n, m, fu)
% 生成 右端向量, n 为大三对角阵阶数, m 子阵阶数, fb(i,j) 为差分方程右端函数
% fu(i,j) 为边界函数
    b = zeros(n*m,1);
    for i = 1:n
        for j = 1:m
            b((i-1)*m+j) = fb(i,j);
        end
        if i == 1 || i == n
            b((i-1)*m+1:i*m) = b((i-1)*m+1:i*m) + fu(i,0);  % 加入 u(i,0)
            b((i-1)*m+1) = b((i-1)*m+1) + fu(0,1);  % 加入 u(0,1)
            b(i*m) = b(i*m) + fu(n,1);  % 加入 u(n,1)
        else
            b((i-1)*m+1) = b((i-1)*m+1) + fu(0,j);  % 加入 u(0,j)
            b(i*m) = b(i*m) + fu(n,j);  % 加入 u(n,j)
        end
    end
end

function [S] = S137_2(i,n)
% 生成 P137(2) 的 主对角线矩阵 S
    global h_137_2; 
    S = zeros(n);
    for j = 1:n
        if j>1 ; S(j,j-1) = -1; end
        S(j,j) = 4 + h_137_2^2*g137_2(i*h_137_2, j*h_137_2);
        if j<n ; S(j,j+1) = -1; end
    end
end

function [X] = p137_2(N)
    global h_137_2; 
    tic
    h_137_2 = 1/N;
    B = -eye(N-1);
    A = createA(@S137_2, B, N-1);
    b = createb(@b137_2, N-1, N-1, @fu137_2);
    X = Solve_GaussSeidel(A,b,rand(1,(N-1)^2),0.00000001);
    toc
end

function [r]=f159_1(x,y)
    r = sin(x*y);
end

function [b] = b159_1(i,j)
% 计算 P159(1) 差分方程组右端向量
    global h_159_1; 
    b = h_159_1^2 * f159_1(i*h_159_1, j*h_159_1)/4;
end

function [u] =fu159_1(i,j)
    global h_159_1; 
    u = (i*h_159_1)^2 + (j*h_159_1)^2;
end

function [S] = S159_1(i,n)
% 生成 P159(1) 的 主对角线矩阵 S
    global h_159_1; 
    S = zeros(n);
    for j = 1:n
        if j>1 ; S(j,j-1) = -1/4; end
        S(j,j) = 1 + h_159_1^2/4;
        if j<n ; S(j,j+1) = -1/4; end
    end
end


function [X] = p159_1(N)
    global h_159_1; 
    tic
    h_159_1 = 1/N;
    B = -eye(N-1)./4;
    A = createA(@S159_1, B, N-1);
    b = createb(@b159_1, N-1, N-1, @fu159_1);
    X = Solve_ConjugateGradient(A,b,rand(1,(N-1)^2),0.0001);
    toc
end

function [X] = p159_1_0(N)
    global h_159_1; 
    tic
    h_159_1 = 1/N;
    B = -eye(N-1)./4;
    A = createA(@S159_1, B, N-1);
    b = createb(@b159_1, N-1, N-1, @fu159_1);
    X = Solve_SOR(A,b,rand(1,(N-1)^2),1.7,0.0001);
    toc
end


function [e]=p160_2(N)
    A = createHilbert(N);
    b = zeros(N,1);
    for i = 1:N; b(i) = sum(A(i,:))/3; end
    X = Solve_ConjugateGradient(A,b,rand(1,N),0.0000001);
    e=norm(X-1/3);
end


function p160_3()
    A = [10 1 2 3 4; 1 9 -1 2 -3; 2 -1 7 3 -5; 3 2 3 12 -1; 4 -3 -5 -1 15];
    b = [12 -27 14 -17 12]';
    x0 = rand(1,5);
    tic
    Solve_Jacobi(A,b,x0,0.00001)'
    toc
    tic
    Solve_GaussSeidel(A,b,x0,0.00001)'
    toc
    tic
    Solve_ConjugateGradient(A,b,x0,0.00001)'
    toc
end