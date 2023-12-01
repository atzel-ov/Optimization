function grad = g(theta, lambda, y, X)

    grad = (1/length(y)) * X * ( 1./( 1+exp( -X'*theta ) )- y) + lambda*theta;

end