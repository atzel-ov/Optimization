clear
close all 
clc

theta = 0:0.01:2*pi;
x1 = -1.5:0.01:2;
y1 = zeros(size(x1));
x_circle = cos(theta);
y_circle = sin(theta);
plot(x_circle, y_circle, 'c', 'LineWidth', 1.3)
hold on
plot(x1, y1, 'b', LineWidth=1.3)
hold on
plot(y1, x1, 'r', LineWidth=1.3)
hold on
plot(x1,-x1+1, 'g', LineWidth=1.3)
hold on
x = 0:0.01:1;
y = zeros(size(x));
plot(x, y, 'm-.', LineWidth=2)
hold on
plot(y, x, 'm-.', LineWidth=2)
hold on
plot(x, -x+1, 'm-.', LineWidth=2)
hold on 
text(0.25,0.25,'$\mathbf{S}$', 'Interpreter', 'latex','FontSize',14)
axis equal
axis([-1.5 2 -1.5 2])
legend('$\|\mathbf{x}\| \leq 1 $', '$x_1 \geq 0$', '$x_2 \geq 0$', '$x_1+x_2\leq 1$', 'Interpreter', 'latex')
