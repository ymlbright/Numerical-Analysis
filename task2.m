function task2()
    clc
    %p98_6(13,5,30)
    %p98_11(30)
    %p99_1()
     p99_2()
     p99_3()
%     [x,y,z]=p137_1(1);
%     [x,y,z]
%     [x,y,z]=p137_1(0.1);
%     [x,y,z]
%     [x,y,z]=p137_1(0.01);
%     [x,y,z]
%     [x,y,z]=p137_1(0.0001);
%     [x,y,z]

end

function p98_6(i, k, n)
    [c,s] = Solve_Givens(0,1);
    A = eye(n);
    A(i,i) = c;
    A(i,k) = -s;
    A(k,i) = s;
    A(k,k) = c;
    x = zeros(n,1);
    y = x;
    x(i) = 1;
    y(k) = 1;
    norm(A'*A-eye(n))
    norm(A*x-y)
end

function p98_11(n)
    A = eye(n);
    A(1,2:n) = rand(1,n-1);
    A(2:n,1) = rand(n-1,1);
    Q = eye(n);
    B = A;
    for i = 1:n-1
        [D,A(i:n,i:n)] = p98_11_p(A(i:n,i:n));
        H = eye(n);
        H(i:n,i:n) = D;
        Q = H*Q;
    end
    norm(Q'*A - B)
    norm(Q'*Q - eye(n))
end


    

function [Q, A] = p98_11_p(A)
    n = length(A);
    Q = eye(n);
    for i = 2:n
        [c, s] = Solve_Givens(A(1,1),A(i,1));
        D = eye(n);
        D(1,1) = c;
        D(1,i) = s;
        D(i,1) = -s;
        D(i,i) = c;
        A = D*A;
        Q = D*Q;
    end
end

function [A] = createTriDiag(d, e, f, n)
    A = zeros(n);
    for i = 1:n
        A(i,i) = e(i);
        if i<n ; A(i,i+1) = f(i); end
        if i>1 ; A(i,i-1) = d(i-1); end
    end
end


function [A] = createHilbert(n)
    A = zeros(n,n);
    for i = 1:n
        for j  = i:n
            A(i,j) = 1/(i+j-1);
            A(j,i) = A(i,j);
        end
    end
end

function p99_1()
    d = 8*ones(83,1);
    e = 6*ones(84,1);
    f = 1*ones(83,1);
    A = createTriDiag(d, e, f, 84);
    b = 15*ones(84,1);
    b(1) = 7;
    b(84) = 14;
    result = zeros(84,2);
    result(:,1) = Solve_Gauss(A, b);
    result(:,2) = Solve_QR(A, b);
    norm(result(:,1)-1)
    norm(result(:,2)-1)
    
    d = 1*ones(99,1);
    e = 10*ones(100,1);
    f = 1*ones(99,1);
    A = createTriDiag(d, e, f, 100);
    b = rand(100,1);
    x1 = Solve_Cholesky(A, b);
    x2 = Solve_QR(A, b);
    norm(A*x1-b)
    norm(A*x2-b)
    
    b = zeros(40,1);
    A = createHilbert(40);
    for i = 1:40
        b(i) = sum(A(i,:));
    end
    x1 = Solve_Cholesky(A, b);
    x2 = Solve_QR(A, b);
    norm(A*x1-b)
    norm(A*x2-b)
end

function p99_2()
    t = -1:0.25:0.75;
    y = [1 0.8125 0.75 1 1.3125 1.75 2.3125]';
    A = zeros(7,3);
    for i = 1:7
        A(i,:) = [t(i)^2 t(i) 1];
    end
    Solve_QR(A,y)
end

function p99_3()
    A = rand(28,12);
    A(:,1) = ones(28,1);
    b = rand(28,1);
    Solve_QR(A,b)
end

function [x, y, z] = p137_1(yipilo)
    n = 100;
    h = 1/n;
    
    d = (yipilo+h)*ones(n-1,1);
    e = -(2*yipilo + h)*ones(n,1);
    f = yipilo*ones(n-1,1);
    A = createTriDiag(d, e, f, n);
    b = (0.5*h^2)*ones(n,1);
    x = 0:h:1-h;
    y = 0.5*(1-exp(-x./yipilo))/(1-exp(-1/yipilo))+0.5*x;
    X = rand(1,n);
    x = norm(Solve_Jacobi(A,b,X,0.00000001)-y');
    y = norm(Solve_GaussSeidel(A,b,X,0.00000001)-y');
    z = norm(Solve_SOR(A,b,X,1.1,0.00000001)-y');
end