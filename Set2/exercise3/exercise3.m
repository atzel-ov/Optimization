clear
close all
clc

%% Parameter Initialization

n = 300; m = 800;

c = randn(n,1);
b = 5*randn(m,1) + 50;
A = 3*randn(m,n);

%% Minimization of f(x) using CVX

cvx_begin

    variable x(n)
    minimize(f(A, b, c, x))

cvx_end

fprintf("============================================================================================\n\n")

%% Contour & Mesh of the Cost function

if(n == 2)
 
x1 = -50:0.1:50;
x2 = -50:0.1:50;

[X1,X2] = meshgrid(x1, x2);

for i = 1:size(X1)
    for j = 1:size(X1)
        x_t = [X1(i,j); X2(i,j)];
        domain_flag = dom_f(A, b, c, x_t);
        if(domain_flag == true)
            cost_function(i,j) = f(A, b, c, x_t);
        else
            cost_function(i,j) = 1E3;
        end
    end
end


figure
mesh(X1,X2,cost_function)
figure
contour(X1,X2,cost_function)
hold on
end


%% Gradient Descent with Backtracking line-search

x0 = zeros(n,1);
epsilon = 0.001; 
alpha = .1; beta = .7;

[x_opt, p_opt1, back_rec] = gradient_method_backtracking(A, b, c, x0, epsilon, alpha, beta);

if (n == 2)
    fprintf("The backtracking line-search solution is x* = (%f,%f)\n\n", x_opt(1), x_opt(2))
    plot(x_opt(1),x_opt(2),'r*')
end

fprintf("The value of the function is p* = %f\n\n", p_opt1)
fprintf("============================================================================================\n\n")

%% Newton-Raphson Algorithm with Backtracking line-search

x0 = zeros(n,1);
epsilon = 0.001; 
alpha = .1; beta = .7;

[x_opt, p_opt2, newton_rec] = newton_method_backtracking(A, b, c, x0, epsilon, alpha, beta);

if (n == 2)
    fprintf("The Newton Algorithm solution is x* = (%f,%f)\n\n", x_opt(1), x_opt(2))
end

fprintf("The value of the function is p* = %f\n\n", p_opt2)
fprintf("============================================================================================\n\n")

%% Convergence Analysis

figure
semilogy(back_rec(1,:),(back_rec(2,:)-p_opt1), 'b', LineWidth=1.5)
hold on
semilogy(newton_rec(1,:),(newton_rec(2,:)-p_opt2), 'r', LineWidth=1.5)
ylabel('$log(f(\mathbf{x}_k) - p_*)$','Interpreter','latex')
xlabel('$k$','Interpreter','latex')
legend('Gradient Method', 'Newton Method')
axis tight
grid on