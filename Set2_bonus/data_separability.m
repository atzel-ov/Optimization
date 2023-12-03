function flag = data_separability(w, b, y, X)

    for i = 1:length(y)
        if(y(i) == 0)
            y(i) = -1;
        end
    end

    flag = 1;
    for i = 1:length(y)
        if(y(i)*(w'*X(:,i)-b) < 0)
            flag = 0;
            break;
        end
    end

end
