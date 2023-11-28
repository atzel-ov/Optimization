clear
close all
clc

%% Construction of Quadratic Function
n = 200;

A = randn(n,n);
[U,S,V] = svd(A);

L = 100;
l = 1;

z = l + (L - l)*rand(n-2,1);

eig_P = [l; L; z];

Lambda = diag(eig_P);

P = U*Lambda*U';     

q = randn(n,1)+0.0001;

c = 2;

%% CVX solution

cvx_begin

    variables x(n);
    minimize(f(P,q,x));
    subject to 
        x'*x - c^2 <= 0;    

cvx_end
fprintf("============================================================================================\n\n")


%% Gradient Algorithm results
if(n == 2)
contour_f(P,q,n); 
theta = linspace(0, 2*pi, 100);
x_circle = c*cos(theta);
y_circle = c*sin(theta);
plot(x_circle, y_circle, 'm--', 'LineWidth', 1.3)
end

x0 = 4*ones(n,1);
epsilon = 10^-6;

[x_opt1, p_opt1, rec1] = gradient_algorithm(P, q, x0, epsilon, @(x)projection_d(c,x));


if (n == 2)
fprintf("The backtracking line-search solution is x* = (%f,%f)\n\n", x_opt1(1), x_opt1(2)); end
fprintf("The value of the function is p* = %f\n\n", p_opt1)
fprintf("============================================================================================\n\n")


%% Accelerated Gradient Algorithm results
if(n == 2)
contour_f(P,q,n); 
theta = linspace(0, 2*pi, 100);
x_circle = c*cos(theta);
y_circle = c*sin(theta);
plot(x_circle, y_circle, 'm--', 'LineWidth', 1.3)
end

x0 = 4*ones(n,1);
epsilon = 10^-6;
beta = (sqrt(L)-sqrt(l))/(sqrt(L)+sqrt(l));

[x_opt2, p_opt2, rec2] = accelerated_gradient_algorithm(P, q, x0, epsilon, beta, @(x)projection_d(c,x));

if (n == 2)
fprintf("The backtracking line-search solution is x* = (%f,%f)\n\n", x_opt2(1), x_opt2(2)); end
fprintf("The value of the function is p* = %f\n\n", p_opt2)
fprintf("============================================================================================\n\n")

%% Convergence Analysis
figure
semilogy(rec1(1,:),abs((rec1(2,:)-p_opt1)), 'g', LineWidth=1.5)
hold on
semilogy(rec2(1,:),abs((rec2(2,:)-p_opt2)), 'm', LineWidth=1.5)
ylabel('$log(f(\mathbf{x}_k) - p_*)$','Interpreter','latex')
xlabel('$k$','Interpreter','latex')
legend('Gradient', 'Accelerated Gradient')
axis tight
grid on