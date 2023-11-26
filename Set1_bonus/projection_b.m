function y = projection_b(a, b, x)

    y = x;

    if(x > a)
        y = min(x, b);
    end

    if(x < b)
        y = max(x, a);
    end

end
