function [b,residuals]=trimmedLinearRegression(Y,X,seed,tolerance);
%function [b,residuals]=trimmedLinearRegression(Y,X,seed,tolerance);

nVars = size(X,2);
constantPos = find(all(X == 1));

if nargin<3, seed = []; end;
if nargin<4, tolerance = []; end;

if isempty(tolerance)
    tolerance = 0.0001; 
end;

constantPos = find(all(X == 1));
if isempty(constantPos)
    X = [ones(size(X,1),1) X];
    nVars = size(X,2);
elseif length(constantPos)>1
    X(:,constantPos(2:end)) = [];
    constantPos = constantPos(1);
end
if ~(constantPos==1)
    tmp = X;
    X(:,1) = X(:,constantPos);
    X(:,2:nVars-1) = X(:,[1:constantPos-1 constantPos+1:nVars]);
    constantPos = 1;
end

if isempty(seed)
    seed = [mean(Y) - sum(mean(X(:,[2:nVars]))); ones(nVars-1,1)];
end
if size(seed,1) == 1
    seed = seed';
end

b_new = seed;
err = Inf;
step = 1;
while (any((err)>tolerance) & step<1000)
    b = b_new;
    residuals = Y - X*b;
    [tmp,keep]=sort(abs(residuals));
    keep = sort(keep(1:round(0.75*length(tmp))));
    b_new = regress(Y(keep),X(keep,:));
    err = b_new - b;
    step = step+1;
end
