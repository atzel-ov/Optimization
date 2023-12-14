function [y_shuffle, X_shuffle] = shuffle_data(y, X)

    A = [y, X'];
    rows = size(A,1);
    P = randperm(rows);
    B = A(P,:);
    y_shuffle = B(:,1);
    X_shuffle = B(:,2:end)';

end