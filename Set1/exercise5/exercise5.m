clear
close all
clc

% C
x = 0.01:0.01:5;

a1 = 2;
f1 = x.^a1;

a2 = -2;
f2 = x.^a2;

a3 = 0.33;
f3 = x.^a3;

figure
hold on
plot(x,f1,'b',LineWidth=1.2)
plot(x,f2,'r',LineWidth=1.2)
plot(x,f3,'g',LineWidth=1.2)
legend('\alpha >= 1','\alpha <= 0','0 <= \alpha <= 1')
hold off
xlabel('x')
ylabel('x^\alpha')
axis([0 5 0 5])
grid on

% D

x1 = -4:0.01:4;
x2 = -4:0.01:4;
[X1,X2] = meshgrid(x1, x2);

for i = 1:size(X1)
    for j = 1:size(X1)
        x_t = [X1(i,j); X2(i,j)];
        f1(i,j) = sqrt(x_t(1)^2+x_t(2)^2);
        f2(i,j) = x_t(1)^2+x_t(2)^2;
    end
end

figure
subplot(1,2,1)
mesh(X1,X2,f1)
xlabel('x_1')
ylabel('x_2')
zlabel('||x||')
grid on
subplot(1,2,2)
mesh(X1,X2,f2)
grid on
xlabel('x_1')
ylabel('x_2')
zlabel('||x||^2')