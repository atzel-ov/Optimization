function quad = f(P, q, x)
    
    quad = 0.5*x'*P*x + q'*x;

end