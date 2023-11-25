function grad = g(A, b, c, x)

    buffer = zeros(size(c,1),1);

    for i = 1:size(b,1)
        buffer = buffer + (1/(b(i)-A(i,:)*x))*A(i,:)';
    end

    grad = c + buffer;

end