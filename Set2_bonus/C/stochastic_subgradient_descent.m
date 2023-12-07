function [theta, fun_value, rec] = stochastic_subgradient_descent(y, X, lambda, theta0, theta_cvx, Ne, Nb, epsilon, k_max)

    theta = theta0;
    fun_value = f(theta, lambda, y, X);
    k = 0;

    e_rec = norm(theta0-theta_cvx);
    k_rec = k;

    while(norm(theta-theta_cvx)>epsilon)

        for i = 1:Ne

            [y_E, X_E] = shuffle_data(y, X);

            for j = 1:length(y)/Nb

                y_B = y_E((j-1)*Nb+1:j*Nb);
                X_B = X_E(:,(j-1)*Nb+1:j*Nb);

                ksi = 0;
                for idx = 1:length(y_B)
                    ksi = ksi + subgrad(y_B(idx),X_B(:,idx),theta);
                end

                ksi = 1/length(y_B)*ksi;
                
                delta_theta = ksi + lambda*theta;

                gamma = 1/(lambda*(k+1));

                theta = theta - gamma*delta_theta;

            end
        end

        k = k + 1;

        fun_value = f(theta, lambda, y_E, X_E);

        e_rec = [e_rec, norm(theta-theta_cvx)];
        k_rec = [k_rec, k];
        
        if(mod(k,5) == 0)
            fprintf("Epoch: %d | f = %f\n", k, fun_value);
        end

        if(k == k_max)
            break;
        end

    end

    rec = [k_rec; e_rec];
end