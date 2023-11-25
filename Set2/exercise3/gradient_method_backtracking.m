function [x, fun_value, back_rec] = gradient_method_backtracking(A, b, c, x0, epsilon, alpha, beta)

    x = x0;
    grad = g(A, b, c, x);
    fun_value = f(A, b, c, x);
    k = 0;

    f_rec = fun_value;
    k_rec = k;
    while(norm(grad) > epsilon)
        delta_x = -grad;
    
        t = 1;
    
        % Feasibility Backtracking
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
    
        % Backtracking line-search
        while(f(A,b,c,x+t*delta_x) > fun_value + alpha*t*grad'*delta_x)
            t = beta*t;
        end
    
        x = x + t*delta_x;
    
        grad = g(A, b, c, x);
        fun_value = f(A, b, c, x);
    
        k = k + 1;
    
        k_rec = [k_rec, k];
        f_rec = [f_rec, fun_value];
        
        fprintf("Iteration: %d | norm_grad = %f | f = %f\n", k, norm(grad), fun_value);
    
    end

    back_rec = [k_rec; f_rec];
end