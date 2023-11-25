clear
close all
clc

x1 = 0:0.01:7;
x2 = 0:0.01:7;
[X1,X2] = meshgrid(x1, x2);

f = fun(X1,X2);

% A
figure
mesh(X1,X2,f)
xlabel('x_1')
ylabel('x_2')
zlabel('f(x_1,x_2)')
grid on

% B
figure
contour(X1,X2,f,20)
xlabel('x_1')
ylabel('x_2')
grid on

% C

x0 = [1.8; 1.8];

f0 = fun(x0(1),x0(2));
df0 = grad_fun(x0(1),x0(2));
d2f0 = hess_fun(x0(1),x0(2));

for i = 1:size(X1)
    for j = 1:size(X1)
        x_t = [X1(i,j); X2(i,j)];
        f_1(i,j) = f0 + df0'*(x_t-x0);
        f_2(i,j) = f0 + df0'*(x_t-x0) + 0.5*(x_t-x0)'*d2f0*(x_t-x0);
    end
end

figure
mesh(X1,X2,f,"EdgeColor",'b')
hold on
mesh(X1,X2,f_1,"EdgeAlpha",0.3)
hold on
mesh(X1,X2,f_2,"EdgeAlpha",0.3)
grid on