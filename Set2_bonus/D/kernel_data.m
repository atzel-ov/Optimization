function [y, X] = kernel_data(M, n, std, r)

    X = zeros(2,n);
    y = zeros(n,1);
    c = zeros(2,M);

    for m = 1:M
        c(:,m) = 5*randn(2,1) + 5*rand(2,1);
    end

    for i = 1:n

        m = unidrnd(M);

        X(:,i) = c(:,m) + std*randn(2,1);

        if(norm(X(:,i)-c(:,m)) <= r)
            y(i) = 1;
        else
            y(i) = -1;
        end

    end

end