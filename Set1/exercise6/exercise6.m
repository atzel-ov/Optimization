clear
close all
clc

A = randn(2,2);
P = A*A';
q = randn(2,1);
r = randn;
x_opt = -inv(P)*q;
disp(x_opt);

x1 = x_opt(1)-5:0.01:x_opt(1)+5;
x2 = x_opt(2)-5:0.01:x_opt(2)+5;

[X1,X2] = meshgrid(x1, x2);

for i = 1:size(X1)
    for j = 1:size(X1)
        x_t = [X1(i,j); X2(i,j)];
        f(i,j) = 0.5*x_t'*P*x_t + q'*x_t + r;
    end
end

figure
mesh(X1,X2,f)
xlabel('x_1')
ylabel('x_2')
zlabel('f(x_1,x_2)')
grid on

figure
hold on
contour(X1,X2,f,40)
plot(x_opt(1),x_opt(2),'r')
text(x_opt(1),x_opt(2),'$\mathbf{x}_*$','Interpreter','latex','FontSize',10,'Position',[x_opt(1)+0.2,x_opt(2)])
xlabel('x_1','Interpreter','latex')
ylabel('x_2','Interpreter','latex')
hold off
grid on