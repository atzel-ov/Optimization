function [y, X] = generate_data(N, n, std, w, b, str)

    x0 = projection_H(eye(N,N), w', b, 4*randn(N,1));   % A point in the hyperplane (computed via a simple projection)
    
    
    gen = randn(n,1);

    if(strcmp(str,'lr'))
        y = (sign(gen) + 1)/2;
    end

    if(strcmp(str,'svm'))
        y = sign(gen);
    end
    

    X = zeros(N,n);

    for i = 1:n
        if(y(i) == 1)
            X(:,i) = x0 + w - std*randn(N,1);
        else
            X(:,i) = x0 - w - std*randn(N,1);
        end
    end

end