clear
close all
clc

%% Data Generation for Hard-SVM

N = 2;              % Dimensions
n = 100;            % Data number
std = 1;          % Standard deviation

w = 4*randn(N,1); b = 4*randn;

flag = 0;
while(flag == 0)
    [y, X] = generate_data(N, n, std, w, b, 'svm');
    flag = data_separability(w, b, y, X);
    if(flag == 1)
        break;
    end
end

if(N == 2)
plot_data(w, b, y, X, 'default');
end

theta = [b; w];                 % Augmented parameters
X_augm = [-ones(1,n); X];       % Augmented data

for i = 1:n
    X_bar(:,i) = y(i)*X_augm(:,i);
end

%% Hard-SVM w/ CVX

cvx_begin   
    variables w(N) b;
    minimize(0.5*(w'*w));
    subject to
        for i = 1:n
            1-y(i)*(w'*X(:,i)-b) <= eps;
        end
cvx_end
fprintf("============================================================================================\n\n")

w_cvx = w; b_cvx = b;

if(N == 2)
plot_data(w_cvx, b_cvx, y, X, 'svm');
end

%% Homogeneous Hard-SVM w/ CVX

cvx_begin
    variables theta(N+1);
    minimize (0.5*(theta'*theta))
    subject to
        for i = 1:n
            1 - y(i)*theta'*X_augm(:,i) <= eps;
        end
cvx_end
fprintf("============================================================================================\n\n")

theta_cvx1 = theta;

if(N == 2)
plot_data(theta_cvx1(2:3), theta_cvx1(1), y, X, 'svm');
end

%% Dual Homogeneous Hard-SVM w/ CVX

cvx_begin
    variables a(n);
    maximize(-0.5*a'*(X_bar'*X_bar)*a + ones(n,1)'*a)
    subject to
        for i = 1:n
            a(i) >= eps;
        end
cvx_end

fprintf("Some of the lagrange multipliers are %f, %f, %f, %f  \n", a(1), a(10), a(20), a(40));
fprintf("============================================================================================\n\n")

theta_cvx2 = 0;
for i = 1:n
    theta_cvx2 = theta_cvx2 + a(i)*y(i)*X_augm(:,i);
end

if(N == 2)
plot_data(theta_cvx2(2:3), theta_cvx2(1), y, X, 'svm');
end

%% Computing the primal through the Dual

theta_cvx1_appr = 0;
for i = 1:n
    theta_cvx1_appr = theta_cvx1_appr + a(i)*y(i)*X_augm(:,i);
end

fprintf("Primal: theta_opt = (%f, %f, %f) \n", theta_cvx1(1), theta_cvx1(2), theta_cvx1(3))
fprintf("Using Eq. (23) the Dual: theta_opt = (%f, %f, %f) \n", theta_cvx1_appr(1), theta_cvx1_appr(2), theta_cvx1_appr(3))