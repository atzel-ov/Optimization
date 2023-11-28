clear
close all
clc

%% Exercise 1/2
x0 = [-1.5; 1];
r = 1;
y = [1; 1];
x = r/(norm(x0-y))*(x0-y)+y;

figure
quiver(0, 0, x0(1), x0(2), 'r', 'LineWidth', 1, 'MaxHeadSize', 0.2)
hold on 
quiver(0, 0, x(1), x(2), 'g-', 'LineWidth', 1, 'MaxHeadSize', 0.2)
theta = linspace(0, 2*pi, 100);
x_circle = y(1) + r*cos(theta);
y_circle = y(2) + r*sin(theta);
plot(x_circle, y_circle, 'b', 'LineWidth', 1.2)
%axis([-2 2 -1.65 1.65])
axis([-2 2.5 -.5 2.5])
axis equal
xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
legend('$\mathbf{x}_0$','$\mathbf{x}$', 'Interpreter', 'latex')
