function sg = subgrad(y, X, theta)

    if(y*theta'*X < 1)
        sg = -y*X;
    else
        sg = 0;
    end

end