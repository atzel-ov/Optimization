function [x, fun_value_x, rec] = accelerated_gradient_algorithm(P, q, x0, epsilon, beta, projection)

    x = x0;
    y = x0;
    fun_value_x = f(P, q, x); grad_x = g(P, q, x);
    fun_value_y = f(P, q, y); grad_y = g(P, q, y);
    t = 1/max(eig(P));
    k = 0;
    
    x_rec = x0;
    f_rec = fun_value_x;
    k_rec = k;

    while(norm(projection(y-t*grad_y)-x)/norm(projection(y-t*grad_y)) > epsilon)
        delta_y = -grad_y;

        x_ = x;
        x = projection(y+t*delta_y);
        y = x + beta*(x - x_);

        grad_x = g(P, q, x); fun_value_x = f(P, q, x);
        grad_y = g(P, q, y); fun_value_y = f(P, q, y);

        k = k + 1;

        x_rec = [x_rec, x];
        k_rec = [k_rec, k];
        f_rec = [f_rec, fun_value_x];

        fprintf("Iteration: %d | norm_grad = %f | f = %f\n", k, norm(grad_x), fun_value_x);

        if(size(x,1) == 2)
            plot(x_rec(1,:),x_rec(2,:),'-ro',LineWidth=0.85)
            hold on
        end

    end

    rec = [k_rec; f_rec];
end