clear
close all
clc

%% Data Generation

N = 2;      % Dimensions
n = 100;    % Data number
std = 1.5;    % Standard deviation

w = 4*randn(N,1); b = 4*randn;

[y, X] = generate_data(N, n, std, w, b, 'lr');

data_separability(w, b, y, X);

theta = [b; w];                 % Augmented parameters
X_augm = [-ones(1,n); X];       % Augmented data
lambda = 10^-2;                 % Regularization parameter

cost = Jr(theta, lambda, y, X_augm);

theta0 = 8*randn(N+1,1);

epsilon = 0.01; k_max = 500;
alpha = 0.1; beta = 0.7;

%% CVX solution

m = N + 1;
cvx_begin
    
    variable theta(m)
    minimize(-(1/length(y)).*(y'*X_augm'*theta) + (1/length(y))*sum(log(1+exp(X_augm'*theta))) + lambda * 0.5*(theta'*theta))

cvx_end
theta_cvx = theta;
fprintf("============================================================================================\n\n")

%% Gradient Descent (GD) w/ Backtracking

[theta_opt1, J_opt1, rec1] = gradient_descent_backtracking(y, X_augm, lambda, theta0, theta_cvx, epsilon, alpha, beta, k_max);

fprintf("The value of the function is p* = %f\n\n", J_opt1)
fprintf("============================================================================================\n\n")

if(N == 2)
data_separability(theta_opt1(2:3), theta_opt1(1), y, X);
end

%% Accelerated Gradient Descent (AGD) w/ Backtracking

[theta_opt2, J_opt2, rec2] = accelerated_gradient_descent(y, X_augm, lambda, theta0, theta_cvx, epsilon, alpha, beta, k_max);

fprintf("The value of the function is p* = %f\n\n", J_opt2)
fprintf("============================================================================================\n\n")

if(N == 2)
data_separability(theta_opt2(2:3), theta_opt2(1), y, X);
end

%% Stochastic Gradient Descent (SGD) 

gamma = 0.3;                    % Learning Rate
epoch_size = size(rec2, 2);     % Epoch size: number of iters od AGD
batch_size = 10;                % Batch size

[theta_opt3, J_opt3, rec3] = stochastic_gradient_descent(y, X_augm, lambda, theta0, theta_cvx, epoch_size, batch_size, epsilon, gamma, k_max);

fprintf("The value of the function is p* = %f\n\n", J_opt2)
fprintf("============================================================================================\n\n")

if(N == 2)
data_separability(theta_opt2(2:3), theta_opt2(1), y, X);
end

%% Convergence Analysis

figure
hold on
semilogy(rec1(1,:), rec1(2,:), 'm', LineWidth=1.3)
semilogy(rec2(1,:), rec2(2,:), 'g', LineWidth=1.3)
semilogy(rec3(1,:), rec3(2,:), 'c', LineWidth=1.3)
hold off
xlabel('$k$', 'Interpreter', 'latex')
ylabel('$\log \|\theta_k - \theta_{cvx}\|$', 'Interpreter', 'latex')
legend('GD', 'AGD', 'SGD')
grid on
axis tight