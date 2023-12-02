function [] = plot_data(w, b, y, X)

figure
hold on
for i = 1:length(y)
    if(y(i) == 1)
        plot(X(1,i), X(2,i), 'bo')
    else
        plot(X(1,i), X(2,i), 'rs')
    end
end

plot_x = [min(X(:))-3,  max(X(:))+3];
plot_y = -(w(1)/w(2))*plot_x + b/w(2);
plot(plot_x, plot_y,'k', LineWidth=1.1)
axis ([min(X(:)) max(X(:)) min(X(:)) max(X(:))])
xlabel('$x_1$', 'Interpreter','latex')
ylabel('$x_2$', 'Interpreter','latex')
hold off