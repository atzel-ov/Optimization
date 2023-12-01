function [x, fun_value_x, rec] = accelerated_gradient_descent(y, X, lambda, x0, x_cvx, epsilon, alpha, beta, k_max)
    
    x = x0; 
    z = x0;
    fun_value_x = Jr(x, lambda, y, X); grad_x = g(x, lambda, y, X);
    fun_value_z = Jr(z, lambda, y, X); grad_z = g(z, lambda, y, X);
    k = 0;
    t = 1;

    e_rec = norm(x0-x_cvx);
    k_rec = k;
    
    %while(norm(x-x_cvx)/norm(x_cvx)>epsilon)
    while(norm(z-t*grad_z-x )/norm(z-t*grad_z)>epsilon)
        delta_z = -grad_z;

        t = 1;
        while(Jr(z + t*delta_z,lambda,y,X) > Jr(z,lambda,y,X)+alpha*t*grad_z'*delta_z)
            t = beta*t;
        end

        x_ = x;
        x = z + t*delta_z;
        z = x + beta*(x - x_);

        fun_value_x = Jr(x, lambda, y, X); grad_x = g(x, lambda, y, X);
        fun_value_z = Jr(z, lambda, y, X); grad_z = g(z, lambda, y, X);

        k = k + 1;

        e_rec = [e_rec, norm(x-x_cvx)];
        k_rec = [k_rec, k];
        
        if(mod(k,5) == 0)
            fprintf("Iteration: %d | norm_grad = %f | Jr = %f\n", k, norm(grad_x), fun_value_x);
        end

        if (k == k_max)
            break
        end

    end

    rec = [k_rec; e_rec];
end