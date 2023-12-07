function fun_val = f(theta, lambda, y, X)

    buffer = 0;

    for i = 1:length(y)
        buffer = buffer + max(0, 1-y(i)*(theta'*X(:,i)));
    end

    fun_val = (lambda/2)*(theta'*theta) + (1/length(y))*buffer;

end