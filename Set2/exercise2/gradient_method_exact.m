function [x, fun_value, exact_rec] = gradient_method_exact(P, q, x0, epsilon)

    %% Plotting
    x_opt = -inv(P)*q;
    p_opt = f(P,q,x_opt);

    if(size(x0,1) == 2)
        x1 = x_opt(1)-7:0.1:x_opt(1)+7;
        x2 = x_opt(2)-7:0.1:x_opt(2)+7;
        [X1,X2] = meshgrid(x1, x2);    
        for i = 1:size(X1)
            for j = 1:size(X1)
                x_t = [X1(i,j); X2(i,j)];
                quad(i,j) = f(P, q, x_t);
            end
        end   
        
        figure
        contour(X1,X2,quad,20)
        axis tight
        hold on
        xlabel('$x_1$','Interpreter','latex')
        ylabel('$x_2$','Interpreter','latex')
    end
    
    %% Algorithm
    % Initial Conditions
    x = x0;
    grad = g(P, q, x);
    fun_value = f(P, q, x);
    k = 0;
    % Recordings Initialization
    x_rec = x0;
    f_rec = fun_value;
    k_rec = k;
    
    while(norm(grad) > epsilon)
        delta_x = -grad;

        t = (norm(grad)^2)/(grad'*P*grad);

        x = x + t*delta_x;

        grad = g(P, q, x);
        fun_value = f(P, q, x);

        k = k + 1;

        x_rec = [x_rec, x];
        k_rec = [k_rec, k];
        f_rec = [f_rec, fun_value];

        fprintf("Iteration: %d | norm_grad = %f | f = %f\n", k, norm(grad), fun_value);

        if(size(x,1) == 2)
            plot(x_rec(1,:),x_rec(2,:),'-ro',LineWidth=0.85)
            hold on
        end
    end
    hold off    

    exact_rec = [k_rec; f_rec];    
end