function [x, fun_value, rec] = gradient_algorithm(P, q, x0, epsilon, projection)

    x = x0;
    fun_value = f(P, q, x); grad = g(P, q, x);
    t = 1/max(eig(P));
    k = 0;

    x_rec = x0;
    f_rec = fun_value;
    k_rec = k;

    while(norm(projection(x-t*grad)-x)/norm(projection(x-t*grad)) > epsilon)
        delta_x = -grad;
        
        x = projection(x+t*delta_x);

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

    rec = [k_rec; f_rec];
end