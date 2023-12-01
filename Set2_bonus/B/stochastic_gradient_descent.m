function [theta, fun_value, rec] = stochastic_gradient_descent(y, X, lambda, theta0, theta_cvx, Ne, Nb, epsilon, gamma, k_max)
  
    theta = theta0;
    fun_value = Jr(theta, lambda, y, X);
    k = 0;

    e_rec = norm(theta0-theta_cvx);
    k_rec = k;

    while(norm(theta-theta_cvx)>epsilon)

        for i = 1:Ne
            A = [y, X'];
            rows = size(A,1);
            P = randperm(rows);
            B = A(P,:);
            y_E = B(:,1);
            X_E = B(:,2:end)';
            for j = 1:length(y)/Nb
                
                y_B = y_E((j-1)*Nb+1:j*Nb);
                X_B = X_E(:,(j-1)*Nb+1:j*Nb);

                delta_theta = -(1/Nb)*g(theta, lambda, y_B, X_B);
                   
                theta = theta + gamma*delta_theta;

            end
        end

        k = k + 1;

        fun_value = Jr(theta, lambda, y_E, X_E);

        e_rec = [e_rec, norm(theta-theta_cvx)];
        k_rec = [k_rec, k];
        
        if(mod(k,5) == 0)
            fprintf("Epoch: %d | Jr = %f\n", k, fun_value);
        end

        if(k == k_max)
            break;
        end

    end

    rec = [k_rec; e_rec];
end