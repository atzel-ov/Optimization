function sub = l(y, x, theta)
    
    sub = 1-y*theta'*x;
    if(1-y*theta'*x < 0)
        sub = 0;
    end
    
end