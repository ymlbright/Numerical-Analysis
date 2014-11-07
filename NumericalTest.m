function NumericalTest( method ,x0, xn, n)
%数值分析测试函数
%NumericalTest( method ,x0, xn, n)
%   method   测试目标函数
%   x0,xn    测试区间
%   n        测试精度(等分份数)

    fprintf(1,'Test Function: y=x^4/5+C/x (x!=0), y0=0 (10 points at most):\n')
    r = method(@df, x0, 0, xn, n);  %进行数值计算
    h = (xn-x0)/n;
    x = f(x0:h:xn);
    fprintf(1,'\t\tx\t\ty\t\t\tyc\t\te\n')
    if n<=10;disp([r,x',x'-r(:,2)]);else d=ceil(n/9);
        disp([r([1:d:n,n],:),x([1:d:n,n])',x([1:d:n,n])'-r([1:d:n,n],2)]);end
    hold on
    plot(r(:,1),r(:,2),'b')
    plot(x0:0.01:xn,f(x0:0.01:xn),'k')
    disp('Convergence oder test:')
    r2 = method(@df, x0, 1, xn, n*2);
    r3 = method(@df, x0, 1, xn, n*4);
    d01 = norm(r2(2*n+1,2)-r(n+1,2),1);
    d02 = norm(r3(4*n+1,2)-r2(2*n+1,2),1);
    d11 = norm(r2(1:2:2*n+1,2)-r(1:n+1,2),1);
    d12 = norm(r3(1:4:4*n+1,2)-r2(1:2:2*n+1,2),1);
    d21 = norm(r2(1:2:2*n+1,2)-r(1:n+1,2),2);
    d22 = norm(r3(1:4:4*n+1,2)-r2(1:2:2*n+1,2),2);
    d31 = norm(r2(1:2:2*n+1,2)-r(1:n+1,2),3);
    d32 = norm(r3(1:4:4*n+1,2)-r2(1:2:2*n+1,2),3);
    fprintf(1,'\tConvergence order (1 norm): %f\n',log(d11/d12)/log(2));
    fprintf(1,'\tConvergence order (2 norm): %f\n',log(d21/d22)/log(2));
    fprintf(1,'\tConvergence order (inf norm): %f\n',log(d31/d32)/log(2));
    fprintf(1,'\tConvergence order (end point): %f\n',log(d01/d02)/log(2))
end

function [r] = df(x,y)
    r = x^3-y/x;
end

function [y] = f(x)
    C = (-1/5)*x(1)^5;
    y = (1/5)*x.^4+C./x;
end
