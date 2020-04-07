function filtered = meanfilter(X,N);

X = X(:)'; % Make X a row matrix
M = length(X);

if N == 1
    filtered = X;
else
    if mod(N,2) == 0
        N = N+1;
    end
    
    indices = repmat(1:M,[N 1]) + repmat((-floor(N/2):floor(N/2))',[1 M]);
    
    minIndex = indices(1,1);
    maxIndex = indices(N,M);
    
    nX = interp1(1:M,X,minIndex:maxIndex,'linear','extrap');
    filtered = mean(nX(indices-minIndex+1));
end