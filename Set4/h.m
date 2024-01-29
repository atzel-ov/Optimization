function hessian = h(x)

    hessian = diag(1./(x.^2));
    
end