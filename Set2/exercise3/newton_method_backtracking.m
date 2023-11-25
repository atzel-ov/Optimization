function [x, fun_value, newton_rec] = newton_method_backtracking(A, b, c, x0, epsilon, alpha, beta)

    x = x0;
    hessian = h(A, b, c, x);
    grad = g(A, b, c, x);
    fun_value = f(A, b, c, x);
    lambda = grad'*inv(hessian)*grad;
    k = 0;    

    f_rec = fun_value;
    k_rec = k;

    while(lambda > 2*epsilon)
        delta_x = -inv(hessian)*grad;

        t = 1;

        while(true)
            domain_flag = dom_f(A, b, c, x+t*delta_x);
            if(domain_flag == 0)
                fprintf('Out of Bounds, now backtracking\n')
                t = beta*t;
                continue
            end
            if(domain_flag == 1)
                break
            end
        end

        while(f(A,b,c,x+t*delta_x) > fun_value + alpha*t*grad'*delta_x)
            t = beta*t;
        end

        x = x + t*delta_x;
        
        hessian = h(A, b, c, x);
        grad = g(A, b, c, x);
        fun_value = f(A, b, c, x);
        lambda = grad'*inv(hessian)*grad;

        k = k + 1;
    
        k_rec = [k_rec, k];
        f_rec = [f_rec, fun_value];

        fprintf("Iteration: %d | norm_grad = %f | f = %f\n", k, norm(grad), fun_value);
    
    end

    newton_rec = [k_rec; f_rec];
end