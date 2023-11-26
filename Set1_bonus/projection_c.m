function y = projection_c(P, A, b, x)
    
    [p, n] = size(A);

    F = [P, A'; A, zeros(p,p)];

    u = [x; b];
    
    iF = inv(F);

    y = iF(1:n,:)*u;

end