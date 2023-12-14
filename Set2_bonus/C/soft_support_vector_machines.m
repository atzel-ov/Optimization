clear
close all
clc

%% Data Generation for Soft-SVM

N = 2;              % Dimensions
n = 400;            % Data number
std = 1.4;          % Standard deviation

w = 4*randn(N,1); b = 4*randn;

[y, X] = generate_data(N, n, std, w, b, 'svm');

if(N == 2)
plot_data(w, b, y, X, 'default');
end

theta = [b; w];                 % Augmented parameters
X_augm = [-ones(1,n); X];       % Augmented data

theta0 = 8*randn(N+1,1);

epsilon = 0.001;
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

%% Batch Size = 1

batch_size = 1;

[theta_opt1, p_opt1, rec1] = stochastic_subgradient_descent(y, X_augm, lambda, theta0, theta_cvx, epoch_size, batch_size);

fprintf("The value of the function is p* = %f\n\n", p_opt1)
fprintf("============================================================================================\n\n")

%% Batch Size = 10

batch_size = 10;

[theta_opt2, p_opt2, rec2] = stochastic_subgradient_descent(y, X_augm, lambda, theta0, theta_cvx, epoch_size, batch_size);

fprintf("The value of the function is p* = %f\n\n", p_opt2)
fprintf("============================================================================================\n\n")

%% Batch Size = 50

batch_size = 50;

[theta_opt3, p_opt3, rec3] = stochastic_subgradient_descent(y, X_augm, lambda, theta0, theta_cvx, epoch_size, batch_size);

fprintf("The value of the function is p* = %f\n\n", p_opt3)
fprintf("============================================================================================\n\n")

%% Batch Size = 100

batch_size = 100;

[theta_opt4, p_opt4, rec4] = stochastic_subgradient_descent(y, X_augm, lambda, theta0, theta_cvx, epoch_size, batch_size);

fprintf("The value of the function is p* = %f\n\n", p_opt4)
fprintf("============================================================================================\n\n")

%% Plotting SSG boundary for the last batch size

if(N == 2)
plot_data(theta_opt1(2:3), theta_opt1(1), y, X, 'svm');
title('$|\mathcal{B}| = 1$', 'Interpreter', 'latex')
plot_data(theta_opt2(2:3), theta_opt2(1), y, X, 'svm');
title('$|\mathcal{B}| = 10$', 'Interpreter', 'latex')
plot_data(theta_opt3(2:3), theta_opt3(1), y, X, 'svm');
title('$|\mathcal{B}| = 50$', 'Interpreter', 'latex')
plot_data(theta_opt3(2:3), theta_opt3(1), y, X, 'svm');
title('$|\mathcal{B}| = 100$', 'Interpreter', 'latex')
end

%% Convergence Analysis

figure
hold on
semilogy(rec1(1,:), rec1(2,:), 'm.-', LineWidth=1.0)
semilogy(rec2(1,:), rec2(2,:), 'go-', LineWidth=1.0)
semilogy(rec3(1,:), rec3(2,:), 'cs-', LineWidth=1.0)
semilogy(rec4(1,:), rec4(2,:), 'b+-', LineWidth=1.0)
hold off
xlabel('$k$', 'Interpreter', 'latex')
ylabel('$\log \|\theta_k - \theta_{cvx}\|$', 'Interpreter', 'latex')
legend('$|\mathcal{B}| = 1$', '$|\mathcal{B}| = 10$', '$|\mathcal{B}| = 50$', '$|\mathcal{B}| = 100$', 'Interpreter', 'latex')
grid on
ylim([0 rec1(2,1)/1.2])
xlim([0 rec1(1,end)/2.1])