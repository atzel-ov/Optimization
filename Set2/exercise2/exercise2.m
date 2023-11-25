clear
close all
clc

%% Construction of Quadratic Function
n = 2;          % Dimensions

A = randn(n,n);
[U,S,V] = svd(A);

l_max = 10;     % Condition number
l_min = 1;

z = l_min + (l_max - l_min)*rand(n-2,1);

eig_P = [l_min; l_max; z];

L = diag(eig_P);

P = U*L*U';     % Construction via EVD

q = randn(n,1)+0.0001;

%% Closed form Solution
x_opt = -inv(P)*q;
p_opt = f(P, q, x_opt);

if (n == 2)
    fprintf("The closed form solution is x* = (%f,%f)\n", x_opt(1), x_opt(2))
end

fprintf("\nThe value of the function is p* = %f\n\n", p_opt)
fprintf("============================================================================================\n\n")

%% CVX Solution

cvx_begin
    
    variable x(n)
    minimize(f(P, q, x))

cvx_end
fprintf("============================================================================================\n\n")


%% Exact Line-search Solution

x0 = 5*ones(n,1);
epsilon = 0.01;

[x_opt, p_optl, exact_rec] = gradient_method_exact(P, q, x0, epsilon);

if (n == 2)
    fprintf("The exact line-search solution is x* = (%f,%f)\n\n", x_opt(1), x_opt(2))
end

fprintf("The value of the function is p* = %f\n\n", p_optl)
fprintf("============================================================================================\n\n")

%% Backtracking Line-search Solution

x0 = 5*ones(n,1);
epsilon = 0.01;
alpha = .1; beta = .7;

[x_opt, p_optb, back_rec] = gradient_method_backtracking(P, q, x0, epsilon, alpha, beta);

if (n == 2)
    fprintf("The backtracking line-search solution is x* = (%f,%f)\n\n", x_opt(1), x_opt(2))
end

fprintf("The value of the function is p* = %f\n\n", p_optb)
fprintf("============================================================================================\n\n")

%% Convergence Analysis

figure
semilogy(exact_rec(1,:),(exact_rec(2,:)-p_optl), 'b', LineWidth=1.5)
hold on
semilogy(back_rec(1,:),(back_rec(2,:)-p_optb), 'r', LineWidth=1.5)
ylabel('$log(f(\mathbf{x}_k) - p_*)$','Interpreter','latex')
xlabel('$k$','Interpreter','latex')
legend('Exact', 'Backtracking')
axis tight
grid on

%% Convergence analysis results for strongly convex functions

ke = (l_max/l_min)*log10((f(P,q,x0)-p_opt)/epsilon);

fprintf("The maximum number of iterations that guarantees solution within accuracy %1.2f is: %1.2f\n\n", epsilon, ke)