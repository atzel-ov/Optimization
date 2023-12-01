function fun_value = Jr(theta, lambda, y, X)

    fun_value = -(1/length(y)).*(y'*X'*theta) + (1/length(y))*sum(log(1+exp(X'*theta))) + lambda * 0.5*norm(theta)^2;

end