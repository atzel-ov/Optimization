function domain_flag = dom_f(A, b, c, x)

    domain_flag = 1;

    for i = 1:size(b,1)
        if(b(i) - A(i,:)*x <= 0)
            domain_flag = 0;
            break;
        end
    end

end