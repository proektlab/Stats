function [beta,residual,residualError] = regressLinearToSigmoid(X,Y,N);
if nargin<3, N = 2; end;
seed = [1 log(2) 1];

for j=1:50,  % multiple repeats to avoid local minima
   params{j} = nlinfit(X,Y,@link,seed+10*(rand(1,3)-0.5));
   residuals{j} = link(params{j},X,N)-Y;
   err(j) = mean(residuals{1}.^2);
end
[residualError,use] = min(err);
beta = params{use};
residual = residuals{use};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yhat = link(param,X,N);
if nargin<3, N = 1.75; end;
yhat = param(1)*(X.^N)./(exp(param(2)).^N + X.^N) + param(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = score(param,X,Y,N);
if nargin<4, N = 2; end;
   err = mean((link(param,X,N)-Y).^2);
return;
