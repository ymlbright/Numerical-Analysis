function task1()

% d = 8*ones(83,1);
% e = 6*ones(84,1);
% f = 1*ones(83,1);
% A = createTriDiag(d, e, f, 84);
% b = 15*ones(84,1);
% b(1) = 7;
% b(84) = 14;
% result  = zeros(84,3);
% result(:,1) = Solve_Gauss(A, b);
% result(:,2) = Solve_Chasing(d, e, f, b);
% result(:,3) = result(:,1) - result(:,2);
% result

for n = 5:50
    k = Estimate_Matrix_Condition_Modinf(createHilbert(n),10);
    fprintf(1,'cond(Hilbert[%d], inf) =\t%f\n',n,k);
end
disp('Error of SolveGauss(An)')
x = rand(30,1);
for n = 5:30
    A = createAn(n);
    b = A*x(1:n);
    X = Solve_Gauss(A, b);
    r = norm(b - A*X, inf);
    ereal = norm(X-x(1:n), inf);
    bb = norm(b, inf);
    k = Estimate_Matrix_Condition_Modinf(A, 10);
    e = k*r/bb;
    fprintf(1,'[n=%d, cond=%d, e=%.16f]:\t%.16f\n',n ,k, e, ereal);
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

function [A] = createAn(n)
A = zeros(n,n);
for i = 1:n
    for j = 1:n
        if j==n
            A(j,i) = 1;
        elseif i>j
            A(j,i) = -1;
        elseif i<j
            A(j,i) = 0;
        else
            A(j,i) = 1;
        end
    end
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