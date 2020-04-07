function [beta,residual,residualError,betaCI] = regressPhaseResponse(X,theta)
% function [beta,residual,residualError] = regressLinearToCircular(x,theta);
%theta = angle(exp(complex(0,1).*theta));
seed = [0.30 15];

params{50} = [];
residuals = params;
err = zeros(50,1);
for j=1:50,  % multiple repeats to avoid local minima
   params{j} = fminsearch(@score,seed+10*(rand(1,2)-0.5),[],X,theta);
   [params{j},R{j},J{j}] = nlinfit(X,theta,@link,params{j});
   residuals{j} = circularDistance(link(params{j},X),theta);
   err(j) = mean(residuals{1}.^2);
end
[residualError,use] = min(err);
beta = params{use};
residual = residuals{use};
betaCI = nlparci(beta,R{use},J{use});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function angle = link(param,X)
 tau_0 = 0.5*(1+tanh(param(1)));
 %tau_1 = max(0.5*(1+tanh(param(2))),eps);
 angle = -2.*pi.*96/46*param(1)*X(:,2).*tanh(param(2).*X(:,1));
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = score(param,X,Y)
   err = mean(circularDistance(link(param,X),Y).^2);
return;
