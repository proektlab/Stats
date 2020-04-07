function [beta,residual,residualError] = regressLinearToCircular(X,theta);
% function [beta,residual,residualError] = regressLinearToCircular(x,theta);

seed = [mean(mean(theta(find(X == max(X)))) - mean(theta(find(X == min(X))))) -3];

for j=1:50,  % multiple repeats to avoid local minima
   params{j} = fminsearch(@score,seed+10*(rand(1,2)-0.5),[],X,theta);
   params{j} = nlinfit(X,theta,@link,params{j});
   residuals{j} = circularDistance(link(params{j},X),theta);
   err(j) = mean(residuals{1}.^2);
end
[residualError,use] = min(err);
beta = params{use};
residual = residuals{use};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function angle = link(param,X);
 angle = param(1)*tanh(param(2)*X);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = score(param,X,Y);
   err = mean(circularDistance(link(param,X),Y).^2);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhosq = circcorr(X,Y);
 % From Jupp and Mardia, 1980. (Biometrika 67:1, 163-173).
 Cx = cos(X);
 Sx = sin(X);
 Cy = cos(Y);
 Sy = sin(Y);
 
 Rcc = correlation(Cx,Cy);
 Rcs = correlation(Cx,Sy);
 Rsc = correlation(Sx,Cy);
 Rss = correlation(Sx,Sy);
 
 R1 = correlation(Cx,Sx);
 R2 = correlation(Cy,Sy);
 
 
 rhosq = (Rcc.^2+Rcs.^2+Rsc.^2+Rss.^2 + 2*(Rcc*Rss - Rcs*Rsc)*R1*R2 ...
     - 2*(Rcc*Rcs + Rsc*Rss)*R2 - 2*(Rcc*Rsc + Rcs*Rss)*R1) ...
     ./ ((1 - R1.^2)*(1 - R2.^2));
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = correlation(X,Y);
 X = X(:); Y = Y(:);
 r = corrcoef(X,Y);
 r = r(1,2);
return;
