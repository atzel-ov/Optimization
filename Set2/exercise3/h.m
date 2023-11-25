function hessian = h(A, b, c, x)

    buffer = zeros(size(c,1),size(c,1));

    for i = 1:size(b,1)
        buffer = buffer + (1/(b(i)-A(i,:)*x)^2)*(A(i,:)'*A(i,:));
    end

    hessian = buffer;

end