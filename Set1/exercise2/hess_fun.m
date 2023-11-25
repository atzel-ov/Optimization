function hess_f = hess_fun(x1,x2)
    
    hess_f = (2./((1+x1+x2)^3))*ones(2,2);

end