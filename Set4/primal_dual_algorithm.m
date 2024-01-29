function [x, fun_value, rec] = primal_dual_algorithm(A, b, c, x0, e_feas, e, alpha, beta, m, mu)

[p, n] = size(A);

x = x0;

t = 1;

lambda = (1/t)*ones(n,1)./x;
v = A'\(lambda-c);

y = [x; lambda; v];

k = 0;
h = x'*lambda;

x_rec = x0;
f_rec = f(c,x);
k_rec = k;

while(1)

    t = mu*m/h;

    F = [zeros(n,n), -eye(n), A'; diag(lambda), diag(x), zeros(n,p); A, zeros(p,n), zeros(p,p)];

    delta_y = -F\r(A, b, c, t, y);
    delta_x = delta_y(1:n);
    delta_l = delta_y(n+1:2*n);
    delta_v = delta_y(2*n+1:2*n+p);

    buffer_l = 1;
    for i = 1:n
        if(-lambda(i)/delta_l(i)<buffer_l && delta_l(i) < 0)
            buffer_l = -lambda(i)/delta_l(i);
        end
    end

    s_max = min(1,buffer_l);

    s = 0.99*s_max;

    x_new = x + s*delta_x;
    %disp(x_new')
    while(not(all(x_new > -10^-15)))
        s = beta*s;
        x_new = x + s*delta_x;
    end

    lambda_new = lambda + s*delta_l;
    v_new = v + s*delta_v;
    y_new = [x_new; lambda_new; v_new];

    while(norm(r(A, b, c, t, y_new)) > (1-alpha*s)*norm(r(A, b, c, t, y)))
        s = beta*s;
        x_new = x + s*delta_x;
        lambda_new = lambda + s*delta_l;
        v_new = v + s*delta_v;
        y_new = [x_new; lambda_new; v_new];
    end

    x = x + s*delta_x;
    lambda = lambda + s*delta_l;
    v = v + s*delta_v;

    y = [x; lambda; v];

    k = k + 1;
    fun_value = f(c,x);
    
    x_rec = [x_rec, x];
    k_rec = [k_rec, k];
    f_rec = [f_rec, fun_value];
    
    
    fprintf("Iteration: %d | f = %f\n", k, fun_value);
    

    if(n == 2)
        plot(x_rec(1,:),x_rec(2,:),'-ro',LineWidth=0.85)
        hold on
    end

    h = x'*lambda;
    
    residual = r(A, b, c, t, y);
    if((norm(residual(2*n+1:2*n+p))<=e_feas) && norm(residual(1:n))<=e_feas && (h < e))
        break
    end

    %warning("off")

    if(k == 1000) break; end

end

rec = [k_rec; f_rec];