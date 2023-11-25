function grad_f = grad_fun(x1,x2)
    
    grad_f = (-1./((1+x1+x2)^2))*ones(2,1);

end