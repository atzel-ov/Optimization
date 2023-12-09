clear
close all
clc

%% Data Generation for Soft-SVM

N = 2;              % Dimensions
n = 200;            % Data number
std = 1.4;          % Standard deviation

w = 4*randn(N,1); b = 4*randn;

[y, X] = generate_data(N, n, std, w, b, 'svm');

if(N == 2)
plot_data(w, b, y, X, 'default');
end

theta = [b; w];                 % Augmented parameters
X_augm = [-ones(1,n); X];       % Augmented data

theta0 = 8*randn(N+1,1);

epsilon = 0.01; k_max = 500; 
lambda = 10^-1;

%% CVX solution

cvx_begin
    variables theta(N+1) ksi(n);
    minimize (0.5*lambda*(theta'*theta) + 1/n*ones(n,1)'*ksi)
    subject to
        for i = 1:n
            1 - ksi(i) - y(i)*theta'*X_augm(:,i) <= 0;
            -ksi(i) <= 0;
        end
cvx_end
fprintf("============================================================================================\n\n")

theta_cvx = theta;

if(N == 2)
plot_data(theta_cvx(2:3), theta_cvx(1), y, X, 'svm')
end

%% Stochastic Sungradient Descent (SSG)

epoch_size = 100;
batch_size = 50;

[theta_opt, p_opt, rec] = stochastic_subgradient_descent(y, X_augm, lambda, theta0, theta_cvx, epoch_size, batch_size, epsilon, k_max);

fprintf("The value of the function is p* = %f\n\n", p_opt)
fprintf("============================================================================================\n\n")

if(N == 2)
plot_data(theta_opt(2:3), theta_opt(1), y, X, 'svm');
end
%% Convergence Analysis

figure
hold on
semilogy(rec(1,:), rec(2,:), 'm', LineWidth=1.3)
hold off
xlabel('$k$', 'Interpreter', 'latex')
ylabel('$\log \|\theta_k - \theta_{cvx}\|$', 'Interpreter', 'latex')
grid on
axis tight

