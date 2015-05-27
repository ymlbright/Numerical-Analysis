function task4()
% 偏微分方程数值解
    clc
    %P32_9;
    %P67_3;
    %P138_10;
    %P138_12;
    P176_8;
    disp('Press any key to exit.')
    pause
    close all
end

function [D, E, F, B] = create_1_3_10(a, b, h, lambda, alpha, q, f)
%构造常微分方程导数边界普通迭代向量 [ -u''(x)+q(x)u(x) = f(x); -u'(a)+lambda1*u(a) = alpha; u'(x)+lambda2*u(b) = beta; a<x<b; ] 
% a,b    - 函数边界
% h      - 步长
% lambda - [lambda1, lambda2]:边值函数系数
% alpha  - [alpha, beta]:边值右端常量
% q,f    - 偏微分方程已知函数
%返回值
% D      - 次对角线向量
% E      - 对角线向量
% F      - 超对角线向量
% B      - 右端向量
    n = floor((a+b)/h);
    D = -1*ones(n,1);
    E = zeros(n+1,1);
    F = D;
    B = zeros(n+1,1);
    E(1) = 1 + lambda(1)*h;
    B(1) = h*alpha(1);
    powerh = h^2;
    for i = 2:n
        a = a + h;
        E(i) = 2 + powerh*q(a);
        B(i) = powerh*f(a); 
    end
    E(n+1) = 1 + lambda(2)*h;
    B(n+1) = h*alpha(2);
end

function [D, E, F, B] = create_1_3_15(a, b, h, lambda, alpha, q, f)
%构造常微分方程导数边界普通迭代向量 [ -u''(x)+q(x)u(x) = f(x); -u'(a)+lambda1*u(a) = alpha; u'(x)+lambda2*u(b) = beta; a<x<b; ] 
% a,b    - 函数边界
% h      - 步长
% lambda - [lambda1, lambda2]:边值函数系数
% alpha  - [alpha, beta]:边值右端常量
% q,f    - 偏微分方程已知函数
%返回值
% D      - 次对角线向量
% E      - 对角线向量
% F      - 超对角线向量
% B      - 右端向量
    n = floor((a+b)/h);
    D = -1*ones(n,1);
    E = zeros(n+1,1);
    F = D;
    B = zeros(n+1,1);
    powerh = h^2;
    E(1) = 1 + lambda(1)*h + powerh/2*q(a);
    B(1) = h*alpha(1) + powerh/2*f(a);
    for i = 2:n
        a = a + h;
        E(i) = 2 + powerh*q(a);
        B(i) = powerh*f(a); 
    end
    a = a + h;
    E(n+1) = 1 + lambda(2)*h + powerh/2*q(a);
    B(n+1) = h*alpha(2) + powerh/2*f(a);
end

function [y] = f_P32_9(x); y = exp(x)*sin(x); end
function [y] = q_P32_9(x); y = 1+sin(x); end
function [y] = u_P32_9(x); y = exp(x); end

function P32_9()
    result = zeros(6,6);
    disp('--------------------------1.3.10 Start-------------------------') 
    U = u_P32_9(0:0.25:1);
    result(6,2:4) = U(2:4);
    h = 0.25;
    for i = 1:5
        result(i,1) = h/2^(i-1);
        [d, e, f, b] = create_1_3_10(0, 1, result(i,1), [1 2], [0 3*exp(1)], @q_P32_9, @f_P32_9);
        u = Solve_Chasing(d, e, f, b);
        result(i,2:4) = [u(1+2^(i-1)) u(1+2*2^(i-1)) u(1+3*2^(i-1))];
        result(i,5) = norm(result(i,2:4)-result(6,2:4), Inf);
        if i>1; result(i,6) = result(i-1,5) / result(i,5); else
            figure(1)
            hold on
            title 'P32函数曲线(绿-1.3.10 蓝-1.3.15 黑-精确值) h=0.25'
            plot(0:0.25:1,U,'k')
            plot(0:0.25:1,u,'g')
            figure(2)
            hold on
            title 'P32误差曲线(绿-1.3.10 蓝-1.3.15)  h=0.25'
            plot(0:0.25:1,abs(U-u'),'g-')
        end
    end
    disp(vpa(result,4))   % 表1.13 填入部分
    disp('--------------------------1.3.10 End---------------------------')
    disp(' ')
    disp('--------------------------1.3.15 Start-------------------------') 
    for i = 1:5
        result(i,1) = h/2^(i-1);
        [d, e, f, b] = create_1_3_15(0, 1, result(i,1), [1 2], [0 3*exp(1)], @q_P32_9, @f_P32_9);
        u = Solve_Chasing(d, e, f, b);
        result(i,2:4) = [u(1+2^(i-1)) u(1+2*2^(i-1)) u(1+3*2^(i-1))];
        result(i,5) = norm(result(i,2:4)-result(6,2:4), Inf);
        if i>1; result(i,6) = result(i-1,5) / result(i,5); else
            figure(1)
            hold on
            plot(0:0.25:1,u,'b')
            figure(2)
            hold on
            plot(0:0.25:1,abs(U-u'),'b-')
        end
    end
    disp(vpa(result,4))   % 表1.13 填入部分
    disp('--------------------------1.3.15 End---------------------------')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = tridiag(a, b, c, n)
    T = b*diag(ones(n,1)) + c*diag(ones(n-1,1),1) + a*diag(ones(n-1,1),-1);
end

function [D, E, F, B] = create_2_2_6(ax, bx, ay, by, hx, hy, f, g)
%构造椭圆型方程Dirichlet边界5点差分迭代向量 [ - ( u[x]''+u[y]'' ) = f(x,y); 边界:u(x,y) = g(x,y); ] 
% ax, bx - x左右边界
% ay, by - y左右边界
% hx     - x方向步长
% hy     - y方向步长
% f      - 右端函数
% g      - 边界函数
%返回值
% D      - 次对角线矩阵
% E      - 对角线矩阵
% F      - 超对角线矩阵
% B      - 右端向量
    n = floor((bx-ax)/hx)-1;
    m = floor((by-ay)/hy)-1;
    E = tridiag(-1/(hx^2), 2*(1/hx^2+1/hy^2), -1/(hx^2), n);
    D = -1/hy^2 * diag(ones(n,1));
    F = D;
    B = zeros(m*n,1);
    mapx = (ax+hx):hx:(bx-hx);
    py = ay;
    for j = 0:m-1;
        py = py + hy;
        B(j*n+1:j*n+n) = f(mapx, py)';
        B(j*n+1) = B(j*n+1) + g(ax, py)/hx^2;
        B(j*n+n) = B(j*n+n) + g(bx, py)/hx^2;
        if j == 0
            B(j*n+1:j*n+n) = B(j*n+1:j*n+n) - D * (g(mapx, ay))';
        elseif j==m-1
            B(j*n+1:j*n+n) = B(j*n+1:j*n+n) - D * (g(mapx, by))';
        end
    end
end

function M = create_map_u(ax, bx, ay, by, hx, hy, f)
% 将网格函数值转换为行优先的列向量
    n = floor((bx-ax)/hx)-1;
    m = floor((by-ay)/hy)-1;
    mapx = (ax+hx):hx:(bx-hx);
    M = zeros(m*n,1);
    for i = 0:m-1
        ay = ay + hy;
        M(i*n+1:i*n+n) = f(mapx, ay)';
    end
end

function draw_3d(ax, bx, ay, by, hx, hy, f, fig, t)
    figure(fig)
    title(t)
    hold on
    n = floor((bx-ax)/hx)-1;
    m = floor((by-ay)/hy)-1;
    Z = zeros(m,n);
    for i = 0:m-1
        Z(i*n+1:i*n+n) =  f(i*n+1:i*n+n);
    end
    surf((ax+hx):hx:(bx-hx),(ay+hy):hy:(by-hy),Z)
    grid on
    view(-37.5,30)
end

function r = f_P67_3(x,y); r=0; end
function r = g_P67_3(x,y); r=(sin(y)+cos(y))*exp(x); end
function r = u_P67_3(x,y); r=(sin(y)+cos(y))*exp(x); end

function P67_3()
    result = zeros(6,9);
    disp('-----------------------------------------2.11 Start----------------------------------------')
    U = create_map_u(0, 1, 0, 1, 0.25, 0.25, @u_P67_3);
    result(6,2:7) = [U(1:3)' U(7:9)'];
    h = 0.25;
    for i = 1:5
        result(i,1) = h/2^(i-1);
        n = floor(1/result(i,1))-1;
        [d, e, f, b] = create_2_2_6(0, 1, 0, 1, result(i,1), result(i,1), @f_P67_3, @g_P67_3);
        u = Solve_ChasingM(d, e, f, n, b);
        basey = n*(n-2^(i-1));
        basex = n*(2^(i-1)-1);
        result(i,2:7) = [u(basex+2^(i-1)) u(basex+2*2^(i-1)) u(basex+3*2^(i-1)) u(basey+2^(i-1)) u(basey+2*2^(i-1)) u(basey+3*2^(i-1))];
        result(i,8) = norm(result(i,2:7)-result(6,2:7), Inf);
        if i>1; result(i,9) = result(i-1,8) / result(i,8); end
        if i == 5
            draw_3d(0, 1, 0, 1, result(i,1), result(i,1), u, 3, 'P67函数数值解曲面 h=1/64');
            U = create_map_u(0, 1, 0, 1, result(i,1), result(i,1), @u_P67_3);
            draw_3d(0, 1, 0, 1, result(i,1), result(i,1), U, 4, 'P67函数精确解曲面 h=1/64');
            draw_3d(0, 1, 0, 1, result(i,1), result(i,1), abs(U-u), 5, 'P67误差曲面 h=1/64');
        end
    end 
    disp(vpa(result,4))   % 表2.11 填入部分
    disp('-----------------------------------------2.11 Start----------------------------------------') 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function U = calc_Z(ax, bx, at, bt, hx, ht, f)
%计算 f 的网格值矩阵
    n = floor((bx-ax)/hx)-1;
    m = floor((bt-at)/ht);
    U = zeros(m,n);
    for i = 1:m
        U(i,:) =  f((ax+hx):hx:(bx-hx),at+i*ht);
    end
end

function U = solve_3_3_5(ax, bx, at, bt, hx, ht, a, f, g, alpha, beta)
%向后Euler法求解抛物型方程Dirichlet边界问题 [ u[t]' - a*u[x]'' = f(x,t); 边界: u(x,at) = g(x); u(ax,t) = alpha(t); u(bx,t) = beta(t); ] 
% ax, bx            - x左右边界
% at, bt            - t左右边界
% hx                - x方向步长
% ht                - t方向步长
% f                 - 右端函数
% f,g,alpha,beta    - 边界函数
%返回值
% U                 - 数值解矩阵，不包含边界
    r = a*ht/hx^2;
    n = floor((bx-ax)/hx)-1;
    m = floor((bt-at)/ht);
    U = zeros(m,n);
    u0 = g((ax+hx):hx:(bx-hx))';
    A = tridiag(-r,1+2*r,-r,n);
    for i = 1:m
        b = u0 + ht*f((ax+hx):hx:(bx-hx),at+i*ht)';
        b(1) = b(1) + r*alpha(at+i*ht);
        b(n) = b(n) + r*beta(at+i*ht);
        u0 = A\b;
        U(i,:) = u0';
    end
end

function r = f_P136_9(x,t); r=-( cos(0.5-t)+2*sin(0.5-t) )*exp(x); end
function r = g_P136_9(x); r = sin(0.5)*exp(x); end
function r = alpha_P136_9(t); r = sin(0.5-t); end
function r = beta_P136_9(t); r = exp(1)*sin(0.5-t); end
function r = real_P136_9(x,t); r = sin(0.5-t)*exp(x); end

function P138_10
    disp('--------------------------3.29 3.30 Start-------------------------')
    result = zeros(10,5);
    result(:,1) = calc_Z(0,1,0,1, 0.5, 0.1, @real_P136_9)';
    hx = 1/100; ht = 1/100;
    u = solve_3_3_5(0,1,0,1, hx, ht, 2, @f_P136_9, @g_P136_9, @alpha_P136_9, @beta_P136_9);
    result(:,2) = u(floor((0.1:0.1:1)/ht),0.5/hx)';
    hx = 1/100; ht = 1/400;
    result(:,3) = abs(result(:,2) - result(:,1));
    u = solve_3_3_5(0,1,0,1, hx, ht, 2, @f_P136_9, @g_P136_9, @alpha_P136_9, @beta_P136_9);
    result(:,4) = u(floor((0.1:0.1:1)/ht),0.5/hx)';
    result(:,5) = abs(result(:,4) - result(:,1));
    disp('    精确解   t=1/10000    误差     t=1/40000   误差')
    disp(vpa(result,4))
    disp('--------------------------3.29 3.30 End---------------------------')

    figure(6)
    title('P138-10 误差曲线图 t=1')
    hold on
    hx = 1/40; ht = 1/1600;
    U = calc_Z(0,1,0,1, hx, 1, @real_P136_9);
    u = solve_3_3_5(0,1,0,1, hx, ht, 2, @f_P136_9, @g_P136_9, @alpha_P136_9, @beta_P136_9);
    plot(hx:hx:1-hx, abs(U-u(length(u),:)), 'r-.')
    hx = 1/20; ht = 1/400;
    U = calc_Z(0,1,0,1, hx, 1, @real_P136_9);
    u = solve_3_3_5(0,1,0,1, hx, ht, 2, @f_P136_9, @g_P136_9, @alpha_P136_9, @beta_P136_9);
    plot(hx:hx:1-hx, abs(U-u(length(u),:)), 'b-*')
    hx = 1/10; ht = 1/100;
    U = calc_Z(0,1,0,1, hx, 1, @real_P136_9);
    u = solve_3_3_5(0,1,0,1, hx, ht, 2, @f_P136_9, @g_P136_9, @alpha_P136_9, @beta_P136_9);
    plot(hx:hx:1-hx, abs(U-u(length(u),:)), 'g-o')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function U = solve_3_5_4(ax, bx, at, bt, hx, ht, a, f, g, alpha, beta)
%Crank-Bicolson求解抛物型方程Dirichlet边界问题 [ u[t]' - a*u[x]'' = f(x,t); 边界: u(x,at) = g(x); u(ax,t) = alpha(t); u(bx,t) = beta(t); ] 
% ax, bx            - x左右边界
% at, bt            - t左右边界
% hx                - x方向步长
% ht                - t方向步长
% f                 - 右端函数
% f,g,alpha,beta    - 边界函数
%返回值
% U                 - 数值解矩阵，不包含边界
    r = a*ht/hx^2;
    n = floor((bx-ax)/hx)-1;
    m = floor((bt-at)/ht);
    U = zeros(m,n);
    u0 = g((ax+hx):hx:(bx-hx))';
    A = tridiag(-r/2,1+r,-r/2,n);
    B = tridiag(r/2,1-r,r/2,n);
    for i = 1:m
        b = B*u0 + ht*f((ax+hx):hx:(bx-hx),at+(i-0.5)*ht)';
        b(1) = b(1) + r*(alpha(at+i*ht)+alpha(at+(i-1)*ht))/2;
        b(n) = b(n) + r*(beta(at+i*ht)+beta(at+(i-1)*ht))/2;
        u0 = A\b;
        U(i,:) = u0';
    end
end

function P138_12
    disp('--------------------------3.33 3.34 3.35 Start----------------------------')
    result = zeros(10,7);
    result(:,1) = calc_Z(0,1,0,1, 0.5, 0.1, @real_P136_9)';
    hx = 1/100; ht = 1/100;
    u = solve_3_5_4(0,1,0,1, hx, ht, 2, @f_P136_9, @g_P136_9, @alpha_P136_9, @beta_P136_9);
    result(:,2) = u(floor((0.1:0.1:1)/ht),0.5/hx)';
    result(:,3) = abs(result(:,2) - result(:,1));
    hx = 1/200; ht = 1/200;
    u = solve_3_5_4(0,1,0,1, hx, ht, 2, @f_P136_9, @g_P136_9, @alpha_P136_9, @beta_P136_9);
    result(:,4) = u(floor((0.1:0.1:1)/ht),0.5/hx)';
    result(:,5) = abs(result(:,4) - result(:,1));
    hx = 1/200; ht = 1/2000;
    u = solve_3_5_4(0,1,0,1, hx, ht, 2, @f_P136_9, @g_P136_9, @alpha_P136_9, @beta_P136_9);
    result(:,6) = u(floor((0.1:0.1:1)/ht),0.5/hx)';
    result(:,7) = abs(result(:,6) - result(:,1));
    disp('   精确解     t=1/100    误差     t=1/200     误差    t=1/20000    误差')
    disp(vpa(result,4))
    disp('--------------------------3.33 3.34 3.35 End------------------------------')
    
    disp('--------------------------3.36 Start----------------------------')
    result = zeros(5,4);
    hx = 0.1; ht = 0.1;
    for i = 1:5
        result(i,1) = hx / 2^(i-1);
        result(i,2) = ht / 2^(i-1);
        U = calc_Z(0,1,0,1, result(i,1), result(i,2), @real_P136_9);
        u = solve_3_5_4(0,1,0,1, result(i,1), result(i,2), 2, @f_P136_9, @g_P136_9, @alpha_P136_9, @beta_P136_9);
        diff = U-u;
        result(i,3) = norm(diff(length(diff),:), Inf); % norm(diff, Inf);
        if i > 1; result(i,4) = result(i-1,3) / result(i,3); end
        if i < 4
            figure(7)
            title('P138-12 误差曲线图 t=1')
            hold on
            plot(result(i,1):result(i,1):1-result(i,1), abs(diff(length(diff),:)), 'b-*')
        end
    end
    disp('      h       t          e       n')
    disp(vpa(result,4))
    disp('--------------------------3.36 End------------------------------')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U] = solve_4_3_10(ax, bx, at, bt, hx, ht, a, f, fa, fb, alpha, beta)
%构造双曲方程隐式差分迭代向量 [ u[t]''-a^2*u[x]'' = f(x,t); 边界: u(x,at) = fa(x); u[t]'(x,at) = fb(x); u(ax,t) = alpha(t); u(ax,t) = alpha(t); u(bx,t) = beta(t);] 
% ax, bx            - x左右边界
% at, bt            - t左右边界
% hx                - x方向步长
% ht                - t方向步长
% f                 - 右端函数
% fa, fb            - t边界函数
% alpha, beta       - x边界函数
%返回值
% U                 - 数值解矩阵，不包含边界
    n = floor((bx-ax)/hx)-1;
    m = floor((bt-at)/ht)-1;
    s = a*ht/hx;
    E = ones(n,1)*(1+s^2);
    D = ones(n-1,1)*(-s^2/2);
    F = D;
    U = zeros(m,n);
    mapx = (ax+hx):hx:(bx-hx);
    U(1,:) = fa(mapx);
    
    
    B = tridiag(s^2/2, 2-s^2, s^2/2, n)*U(1,:)';
    B = B + ht^2*f(mapx,at) + 2*ht*fb(mapx)';
    B(1) = B(1) + s^2/2*(alpha(at+ht) + alpha(at));
    B(n) = B(n) + s^2/2*(beta(at+ht) + beta(at));
    U(2,:) = Solve_Chasing(D, E+1, F, B)';
    
    for j = 2:m+1;
        at = at + ht;
        B = tridiag(s^2/2, -1-s^2, s^2/2, n)*U(j-1,:)';
        B = B + ht^2*f(mapx,at) + 2*U(j,:)';
        B(1) = B(1) + s^2/2*(alpha(at+ht) + alpha(at-ht));
        B(n) = B(n) + s^2/2*(beta(at+ht) + beta(at-ht));
        U(j+1,:) = Solve_Chasing(D, E, F, B)';
    end
    U = U(2:m+2,:);
end

function r = f_P175_7(x,t); r = ((t.^2-x.^2).*sin(x.*t))'; end
function r = fa_P175_7(x); r = zeros(length(x),1); end
function r = fb_P175_7(x); r = x; end
function r = alpha_P175_7(t); r = zeros(length(t),1); end
function r = beta_P175_7(t); r = sin(t); end
function r = real_P175_7(x,t); r = sin(x.*t); end

function P176_8
    U = calc_Z(0,1, 0,1, 0.1, 0.1, @real_P175_7);
    u = solve_4_3_10(0,1,0,1,0.1,0.1,1, @f_P175_7, @fa_P175_7, @fb_P175_7, @alpha_P175_7, @beta_P175_7);
    draw_3d(0,1,0,1, 0.1, 0.1, U, 8, 'P176.10 数值解与精确解图像')
    draw_3d(0,1,0,1, 0.1, 0.1, u, 8, 'P176.10 数值解与精确解图像')
    disp('--------------------------4.14 4.15 Start-------------------------')
    result = zeros(10,5);
    result(:,1) = calc_Z(0,1,0,1, 0.5, 0.1, @real_P175_7)';
    hx = 1/100; ht = 1/200;
    u = solve_4_3_10(0,1,0,1,hx,ht,1, @f_P175_7, @fa_P175_7, @fb_P175_7, @alpha_P175_7, @beta_P175_7);
    result(:,2) = u(floor((0.1:0.1:1)/ht),0.5/hx)';
    result(:,3) = abs(result(:,2) - result(:,1));
    hx = 1/200; ht = 1/100;
    u = solve_4_3_10(0,1,0,1,hx,ht,1, @f_P175_7, @fa_P175_7, @fb_P175_7, @alpha_P175_7, @beta_P175_7);
    result(:,4) = u(floor((0.1:0.1:1)/ht),0.5/hx)';
    result(:,5) = abs(result(:,4) - result(:,1));
    disp('  精确解     s=0.5    误差       s=2      误差')
    disp(vpa(result,4))
    disp('--------------------------4.14 4.15 End---------------------------')
    
    disp('--------------------------4.16 Start----------------------------')
    result = zeros(6,4);
    hx = 0.1; ht = 0.1;
    for i = 1:6
        result(i,1) = hx / 2^(i-1);
        result(i,2) = ht / 2^(i-1);
        U = calc_Z(0,1, 0,1, result(i,1), result(i,2), @real_P175_7);
        u = solve_4_3_10(0,1,0,1,result(i,1),result(i,2),1, @f_P175_7, @fa_P175_7, @fb_P175_7, @alpha_P175_7, @beta_P175_7);
        diff = U-u;
        result(i,3) = norm(diff(length(diff),:), Inf); % norm(diff, Inf);
        if i > 1; result(i,4) = result(i-1,3) / result(i,3); end
    end
    disp('      h         t          e       n')
    disp(vpa(result,4))
    disp('--------------------------4.16 End------------------------------')
    
    h = 1/10;
    for i = 1:3
        U = calc_Z(0,1, 0,1, h, h, @real_P175_7);
        u = solve_4_3_10(0,1,0,1,h,h,1, @f_P175_7, @fa_P175_7, @fb_P175_7, @alpha_P175_7, @beta_P175_7);
        draw_3d(0,1,0,1, h, h, abs(U-u), 9, 'P176.10 0.1,0.05,0.025 三种步长下的误差曲面')
        h = h/2;
    end 
end
