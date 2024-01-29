function [x, fun_value, rec] = interior_point_method(A, b, c_, x0, e_1, e_2, alpha, beta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M-file that implements an interior-point algorithm for the           %
% minimization of a linear function with a linear inequality           % 
% constraint and positivity constraints                                %
%                                                                      %
% Objective function: c^T x                                            %
%                                                                      %
% Inequality constraint:  a_1^T x - b_1 \le 0                          %
%                         x_i >= 0, i=1, 2                             %
%                                                                      %
% Method of solution: Interior Point                                   %
%                                                                      %
% A. P. Liavas, Jan. 24, 2012                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialization 
threshold_1 = e_1;              % should be small !!!
threshold_2 = e_2;              % should be small !!!
h = 10^(-0);                    % "small" parameter for numerical computation of grad, Hess
mu = 2;                         % parameter for increasing t  

% Parameter of linear Cost Function f(x) = c^T x 
c = c_;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Interior point method                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x_init(:,1) = x0;                           % start from a feasible point
t = 1;                                      % initialize parameter t of interior point methods
fun_value = f(c,x_init(:,1));
f_rec = fun_value;
outer_iter = 1;
X(:,outer_iter) = x_init;
k_rec = outer_iter;
x_rec = x0;

l_x = 1000;

while (1)                                     % Outer loop 
     
     x(:,1) = X(:,outer_iter);                 % Start new optimization from the previous solution ...
     inner_iter = 1; 
     while ( 1 )                              % Inner loop 
         
         grad = h* (g(t, c, x(:,inner_iter)));                           % gradient at x
         Hess = h*diag(1./(x(:,inner_iter).^2));                           % Hessian at x
         w = -inv(A*inv(Hess)*A')*A*h*inv(Hess)*grad;
         Dx_Nt = -inv(Hess)*(grad+A'*w);% Newton step
         l_x_ = l_x;
         l_x = sqrt(Dx_Nt'*Hess*Dx_Nt);
         %disp(l_x)% Newton decrement
         if (l_x^2/2 <= threshold_1 || abs(l_x - l_x_)<=10^-7) break; end    % Newton iterations termination condition
         
         % !!! Check feasibility !!!
         tau = 1;
         x_new = x(:,inner_iter) + tau * Dx_Nt;
         while ( not(all(x_new > 0)) )        
             tau = beta * tau; 
             x_new = x(:,inner_iter) + tau * Dx_Nt;
         end 
         % !!!! x_new is FEASIBLE here !!!

         % Backtracking
         while ( t*f(c, x_new) + phi(x_new) > t*f(c, x(:,inner_iter)) + phi(x(:,inner_iter)) + alpha*tau*grad'*Dx_Nt)
              tau = beta * tau; 
              x_new = x(:,inner_iter) + tau * Dx_Nt;
         end

         % Update x
         x(:,inner_iter+1) = x(:,inner_iter) + tau * Dx_Nt;    % Update x
         inner_iter = inner_iter + 1;
     end

     X(:,outer_iter+1) = x(:,inner_iter);           % X: solution of optimization problem      
     if (1/t < threshold_2) break; end              % Algorithm termination condition
     

     fprintf("Iteration: %d | f = %f\n", outer_iter, fun_value);
     warning("off")

     if(size(x(:,inner_iter),1) == 2)
        plot(x_rec(1,:),x_rec(2,:),'-ro',LineWidth=0.85)
        hold on
     end

     outer_iter = outer_iter + 1;
     t = t * mu;

     fun_value = f(c,x(:,inner_iter));
        
     x_rec = [x_rec, x(:,inner_iter)];
     k_rec = [k_rec, outer_iter];
     f_rec = [f_rec, fun_value];

     
end

rec = [k_rec; f_rec];

figure
plot(X(1,:), X(2,:))