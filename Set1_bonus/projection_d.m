function y = projection_d(c, x)

    y = x;

    if(norm(x) > c)
        y = (c/norm(x))*x;
    end

end

