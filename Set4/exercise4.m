clear
close all
clc

%% Parameter Initialization

n = 2;
p = 1;

A = 2*randn(p,n) + 2;
x0 = projection(1.0*randn(n,1)+3*ones(n,1))+ones(n,1);
b = A*x0;

c = 5*randn(n,1) + 1;

alpha = 0.1; beta = 0.7;    % Backtracking parameters

%% CVX Solution

cvx_begin

    variables x(n);
    minimize(f(c,x));
    subject to
        A*x == b;
        x >= eps;

cvx_end

x_star = x;
p_star = f(c,x_star);

nonzero_elements = 0;
for i = 1:length(x_star)
    if(x_star(i) > 10^-5)
        nonzero_elements = nonzero_elements+1;
    end
end


%% Feasibility Checking - Phase 1

cvx_begin quiet

    variables x(n) s;
    minimize(s);
    subject to 
        -x <= s;
        A*x == b;

cvx_end

if(s >= -0.01 || cvx_optval == -Inf)
    fprintf("------------------------------------------\n")
    fprintf("| The problem is infeasible or unbounded |\n")
    fprintf("------------------------------------------\n\n\n")
    return
else 
    fprintf("---------------------------\n")
    fprintf("| The problem is feasible |\n")
    fprintf("---------------------------\n\n\n")
end

%% Interior Point Method

e_1 = 10^-0; e_2 = 10^-7;

if(n==2)
contour_f(c,n,x_star);
x_ = -2:0.1:max(x_star)+5;
y_ = -(A(1)/A(2))*x_ + b/A(2);
plot(x_, y_, 'm', LineWidth=1)
axis([-1 max(x_star)+4 -1 max(x_star)+4])
end

[x_int, fun_value_int, int_rec] = interior_point_method(A, b, c, x0, e_1, e_2, alpha, beta);

fprintf("============================================================================================\n\n")

%% Primal-Dual Newton Algorithm 

mu = 10; m = 1;
e_feas = 10^-5; e = 10^-7;

if(n==2)
contour_f(c,n,x_star);
x_ = -2:0.1:max(x_star)+5;
y_ = -(A(1)/A(2))*x_ + b/A(2);
plot(x_, y_, 'm', LineWidth=1)
axis([-1 max(x_star)+4 -1 max(x_star)+4])
end

[x_pd, fun_value_pd, pd_rec] = primal_dual_algorithm(A, b, c, x0+randn(n,1), e_feas, e, alpha, beta, m, mu);

fprintf("============================================================================================\n\n")
