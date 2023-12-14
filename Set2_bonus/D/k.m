function kernel = k(x, y, str)

    if(strcmp(str, 'Gaussian'))
        l = 1;
        kernel = exp((-norm(x-y).^2)/(2*l^2));
    end

end