clear
close all
clc

%% Data Generation

M = 4;
n = M*100;
std = 0.9;

r = std;

[y, X] = kernel_data(M, n, std, r);

plot_data(zeros(2,1), 0, y, X, 'default');

X_augm = [-ones(1,n); X];
lambda = 10^-1;

for i = 1:n
    for j = 1:n
        K(i,j) = y(i)*y(j)*k(X_augm(:,i), X_augm(:,j), 'Gaussian');
    end
end

%% Dual Homogeneous soft-SVM solution w/ CVX

cvx_begin

    variables a(n);
    buffer = 0;
    maximize(ones(1,n)*a - (0.5/lambda)*a'*K*a);
    subject to 
        a >= eps;
        a <= ones(n,1)/n;

cvx_end

a_opt = a;

%% Grid Generation

x = linspace(min(X(:))-0.3,max(X(:))+0.3,50);

figure
hold on
for i = 1:length(x)
    for j = 1: length(x)
        
        for idx = 1:n
            kernel(idx) = k(X(:,idx),[x(i);x(j)], 'Gaussian');
        end

        y_grid = sign(kernel*(y.*a_opt));
        if(y_grid == 1)
            plot(x(i), x(j), 'ro')
        else 
            plot(x(i),x(j), 'b+')
        end
    end
end