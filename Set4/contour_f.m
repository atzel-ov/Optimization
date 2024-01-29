function [] = contour_f(c,n,x_star)

if(n == 2)
    x1 = -2:0.1:max(x_star)+6;
    x2 = -2:0.1:max(x_star)+6;
    [X1,X2] = meshgrid(x1, x2);
    for i = 1:size(X1)
        for j = 1:size(X1)
            x_t = [X1(i,j); X2(i,j)];
            quad(i,j) = f(c, x_t);
        end
    end   
        
    figure
    contourf(X1,X2,quad,20)
    colormap( flipud(gray(256)) )
    hold on
    plot(x_star(1),x_star(2),'r*')
    axis tight
    xlabel('$x_1$','Interpreter','latex')
    ylabel('$x_2$','Interpreter','latex')
end