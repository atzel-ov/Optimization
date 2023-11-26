function [] = contour_f(P,q,n)

if(n == 2)
    x1 = -7:0.1:+7;
    x2 = -7:0.1:+7;
    [X1,X2] = meshgrid(x1, x2);
    for i = 1:size(X1)
        for j = 1:size(X1)
            x_t = [X1(i,j); X2(i,j)];
            quad(i,j) = f(P, q, x_t);
        end
    end   
        
    figure
    contour(X1,X2,quad,20)
    hold on
    axis tight
    xlabel('$x_1$','Interpreter','latex')
    ylabel('$x_2$','Interpreter','latex')
end
