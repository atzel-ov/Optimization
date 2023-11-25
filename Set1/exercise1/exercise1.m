clear
close all
clc

x = 0:0.01:20;

f = fun(x);

idx = 1;
figure
for x0 = 1:2:8
    
    f0 = fun(x0);
    df0 = fun_1d(x0);
    d2f0 = fun_2d(x0);
    
    f_1 = f0 + df0*(x-x0);
    f_2 = f0 + df0*(x-x0) + 0.5*d2f0*(x-x0).^2;

    subplot(2,2,idx)
    hold on
    plot(x,f,'b',LineWidth=1.2)
    plot(x,f_1,'r',LineWidth=1.2)
    plot(x,f_2,'g',LineWidth=1.2)
    legend('f(x)', 'f_(_1_)(x)', 'f_(_2_)(x)')
    axis([0 10 0 1.5])
    grid on
    title(sprintf('x_0 = %d', x0),"FontSize",8)
    xlabel('x')
    ylabel('f(x)')
    hold off
    
    idx = idx + 1;
end