function [x, fun_value, rec] = gradient_descent_backtracking(y, X, lambda, x0, x_cvx, epsilon, alpha, beta, k_max)
    
    x = x0;
    fun_value = Jr(x, lambda, y, X); grad = g(x, lambda, y, X);
    k = 0;

    e_rec = norm(x0-x_cvx);
    k_rec = k;

    while(norm(grad) > epsilon)
        delta_x = -grad;

        t = 1;
        while(Jr(x + t*delta_x, lambda, y, X) > fun_value + alpha*t*grad'*delta_x)
            t = beta*t;
        end

        x = x + t*delta_x;

        fun_value = Jr(x, lambda, y, X); grad = g(x, lambda, y, X);

        k = k + 1;

        e_rec = [e_rec, norm(x-x_cvx)];
        k_rec = [k_rec, k];

        if(mod(k,5) == 0)
            fprintf("Iteration: %d | norm_grad = %f | Jr = %f\n", k, norm(grad), fun_value);
        end

        if (k == k_max)
            break
        end

    end

    rec = [k_rec; e_rec];
end