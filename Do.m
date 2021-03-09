function r = Do(tau, X)
    % shrinkage operator for singular values
    [U,S,V] = svd(X);
    r = U*So(tau, S)*V';
end
