function [beta,residual,residualError] = regressLinearToArchimedes(X,Y,N);
seed = [1 1 0];

for j=1:50,  % multiple repeats to avoid local minima
   params{j} = nlinfit(X,Y,@link,seed+10*(rand(1,3)-0.5));
   residuals{j} = link(params{j},X)-Y;
   err(j) = mean(residuals{1}.^2);
end
[residualError,use] = min(err);
beta = params{use};
residual = residuals{use};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yhat = link(param,X);
yhat = param(1) + param(2)*log(X) + param(3)*X;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
