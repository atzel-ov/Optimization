clear
close all
clc

%% Data Generation

N = 2;          % Dimensions
n = 100;        % Data number
std = 1.2;      % Standard deviation

w = 4*randn(N,1); b = 4*randn;

[y, X] = generate_data(N, n, std, w, b, 'lr');

flag = data_separability(w, b, y, X);
if (flag == 0)
    fprintf("Data are not separable. \n")
else
    fprintf("Data are separable. \n")
end

if(N == 2)
plot_data(w, b, y, X, 'default');
end

theta = [b; w];                  % Augmented parameters
X_augm = [-ones(1,n); X];        % Augmented data
lambda = 10^-1;                  % Regularization parameter

cost = Jr(theta, lambda, y, X_augm);

theta0 = 8*randn(N+1,1);         % Initial condition

epsilon = 0.001; k_max = 500;    % Stopping condition parameters
alpha = 0.1; beta = 0.7;         % Backtracking parameters

%% CVX solution

m = N + 1;
cvx_begin
    
    variable theta(m)
    minimize(-(1/length(y)).*(y'*X_augm'*theta) + (1/length(y))*sum(log(1+exp(X_augm'*theta))) + lambda * 0.5*(theta'*theta))

cvx_end
theta_cvx = theta;
fprintf("============================================================================================\n\n")

%% CVX solution for various lambdas // Run separably from the rest of the script as cvx is slow; 50 samples are used.
lstr = {'\lambda'};
fprintf("Press 1 to get the results for various lambdas\n")
f = input(' ');
if(f == 1)
m = N+1;
for lambda = [1, 0.1, 0.01, 0.001, 0]
cvx_begin
    
    variable theta(m)
    minimize(-(1/length(y)).*(y'*X_augm'*theta) + (1/length(y))*sum(log(1+exp(X_augm'*theta))) + lambda * 0.5*(theta'*theta))

cvx_end
plot_data(theta(2:3), theta(1), y, X, 'default');
str = strcat(lstr,' = ', num2str(lambda));
title(str)
fprintf("============================================================================================\n\n")
end
end
lambda = 10^-1;  
%% Gradient Descent (GD) w/ Backtracking

[theta_opt1, J_opt1, rec1] = gradient_descent_backtracking(y, X_augm, lambda, theta0, theta_cvx, epsilon, alpha, beta, k_max);

fprintf("The value of the function is p* = %f\n\n", J_opt1)
fprintf("============================================================================================\n\n")

if(N == 2)
plot_data(theta_opt1(2:3), theta_opt1(1), y, X, 'default');
end

%% Accelerated Gradient Descent (AGD) w/ Backtracking

[theta_opt2, J_opt2, rec2] = accelerated_gradient_descent(y, X_augm, lambda, theta0, theta_cvx, epsilon, alpha, beta, k_max);

fprintf("The value of the function is p* = %f\n\n", J_opt2)
fprintf("============================================================================================\n\n")

if(N == 2)
plot_data(theta_opt2(2:3), theta_opt2(1), y, X, 'default');
end

%% Stochastic Gradient Descent (SGD) 

epoch_size = size(rec2, 2);     % Epoch size: number of iters of AGD
batch_size = 10;                % Batch size

[theta_opt3, J_opt3, rec3] = stochastic_gradient_descent(y, X_augm, lambda, theta0, theta_cvx, epoch_size, batch_size);

fprintf("The value of the function is p* = %f\n\n", J_opt3)
fprintf("============================================================================================\n\n")

if(N == 2)
plot_data(theta_opt3(2:3), theta_opt3(1), y, X, 'default');
end

%% Convergence Analysis

figure
semilogy(rec1(1,:), rec1(2,:), 'm', LineWidth=1.1)
hold on
semilogy(rec2(1,:), rec2(2,:), 'g', LineWidth=1.1)
hold on
semilogy(rec3(1,:), rec3(2,:), 'c', LineWidth=1.1)
xlabel('$k$', 'Interpreter', 'latex')
ylabel('$\log \|\theta_k - \theta_{cvx}\|$', 'Interpreter', 'latex')
legend('GD', 'AGD', 'SGD')
grid on
axis tight