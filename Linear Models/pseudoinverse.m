function pinv = pseudoinverse(X)
% function pinv = pseudoinverse(X);
% 
% Returns the pseudoinverse of X.
[U,S,V] = svd(X);
Sinv = zeros(size(S));
Sinv(find(S)) = (1./S(find(S)));
pinv = V*Sinv'*U';