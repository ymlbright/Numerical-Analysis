function [ ret ] = ModEuler( f, x0, y0, xn, n )
%改进欧拉方法数值积分
%[ ret ] = ModEuler( f, x0, y0, xn, n )
%   f  积分函数
%   x0 积分起始点
%   y0 函数初始值
%   xn 积分终止点
%   n  均分份数(精度)
%返回值
%   ret [ x0 y0; x1 y1; ...; xn yn;]

    ret = zeros(n+1,2);
    h = (xn - x0)/n;
    ret(1,2) = y0; 
    ret(:,1) = (x0:h:xn)';
    for i = 2:n+1
        p = ret(i-1,2) + h*f(ret(i-1,1),ret(i-1,2));
        c = ret(i-1,2) + h*f(ret(i,1),p);
        ret(i,2) = (p+c)/2;
    end
end

