function [theta, fun_value, rec] = stochastic_gradient_descent(y, X, lambda, theta0, theta_cvx, Ne, Nb)
  
    theta = theta0;
    fun_value = Jr(theta, lambda, y, X);
    k = 0;

    e_rec = norm(theta0-theta_cvx);
    k_rec = k;
    
    for i = 1:Ne

        [y_E, X_E] = shuffle_data(y, X);

        for j = 1:length(y)/Nb

            y_B = y_E((j-1)*Nb+1:j*Nb);
            X_B = X_E(:,(j-1)*Nb+1:j*Nb);

            delta_theta = -g(theta, lambda, y_B, X_B);

            gamma = 2/(lambda*(k+1));
                   
            theta = theta + gamma*delta_theta;

        end

    k = k + 1;

    fun_value = Jr(theta, lambda, y, X);
    
    e_rec = [e_rec, norm(theta-theta_cvx)];
    k_rec = [k_rec, k];
            
            
    if(mod(k,5) == 0)
        fprintf("Epoch: %d | f = %f\n", k, fun_value);
    end

    end

    rec = [k_rec; e_rec];
end