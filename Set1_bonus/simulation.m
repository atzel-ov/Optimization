clear
close all
clc

%% Construction of Quadratic Function
n = 2;
p = 3;

A = randn(n,n);
[U,S,V] = svd(A);

L = 10;
l = 1;

z = l + (L - l)*rand(n-2,1);

eig_P = [l; L; z];

Lambda = diag(eig_P);

P = U*Lambda*U';     

q = randn(n,1)+0.0001;

a = 1*ones(n,1); b = 2*ones(n,1);
A = 2*rand(p,n); bc = 2*rand(p,1);
c = 2;

%% CVX solution

cvx_begin

    variables x(n);
    minimize(f(P,q,x));
    subject to 
        %-x <= 0;            % Condition a
        %x >= a;             % Condition b
        %x <= b;
        A*x == bc;        % Condition c
        %x'*x - c^2 <= 0;    % Condition d 

cvx_end
fprintf("============================================================================================\n\n")

%% Gradient Algorithm results

%% Projection onto Set (Sa)
if(n == 2)
contour_f(P,q,n); end

x0 = 4*ones(n,1);
epsilon = 10^-6;

[x_opt, p_opta, rec_a] = gradient_algorithm(P, q, x0, epsilon, @projection_a);


if (n == 2)
fprintf("The backtracking line-search solution is x* = (%f,%f)\n\n", x_opt(1), x_opt(2)); end
fprintf("The value of the function is p* = %f\n\n", p_opta)
fprintf("============================================================================================\n\n")


%% Projection onto Set (Sb)
if(n == 2)
contour_f(P,q,n); end

x0 = 4*ones(n,1);
epsilon = 10^-6;

[x_opt, p_optb, rec_b] = gradient_algorithm(P, q, x0, epsilon, @(x)projection_b(a,b,x));


if (n == 2)
fprintf("The backtracking line-search solution is x* = (%f,%f)\n\n", x_opt(1), x_opt(2)); end
fprintf("The value of the function is p* = %f\n\n", p_optb)
fprintf("============================================================================================\n\n")


%% Projection onto Set (Sd)
if(n == 2)
contour_f(P,q,n); end

x0 = 4*ones(n,1);
epsilon = 10^-6;

[x_opt, p_optc, rec_c] = gradient_algorithm(P, q, x0, epsilon, @(x)projection_c(eye(n,n), A, bc, x));


if (n == 2)
fprintf("The backtracking line-search solution is x* = (%f,%f)\n\n", x_opt(1), x_opt(2)); end
fprintf("The value of the function is p* = %f\n\n", p_optc)
fprintf("============================================================================================\n\n")


%% Projection onto Set (Sd)
if(n == 2)
contour_f(P,q,n); end

x0 = 4*ones(n,1);
epsilon = 10^-6;

[x_opt, p_optd, rec_d] = gradient_algorithm(P, q, x0, epsilon, @(x)projection_d(c,x));


if (n == 2)
fprintf("The backtracking line-search solution is x* = (%f,%f)\n\n", x_opt(1), x_opt(2)); end
fprintf("The value of the function is p* = %f\n\n", p_optd)
fprintf("============================================================================================\n\n")