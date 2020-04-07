function [beta,R2] = SHVDregression(x,y,sx,sy);
% function [beta,R2] = SHVDregression(x,y,sx,sy);

if nargin<3, sx = []; end;
if nargin<4, sy = []; end;


correctionFactor = MADn(x(:))./MADn(y(:)); % This should set the slope to one
yprime = y*correctionFactor;

correctedBeta = fminsearch(@gof,[0 1],[],x,yprime,sx,sy);
beta = correctedBeta./correctionFactor;
yhat = beta(1) + beta(2)*x;
residuals = y - yhat;
R2 = 1-var(residuals)./var(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function score = gof(beta,x,y,sx,sy);

[yhat,residuals] = fit(beta,x,y);
if ~isempty(sx), residuals = residuals./max(sx,eps); end;
if ~isempty(sy), residuals = residuals./max(sy,eps); end;

W = (1 + beta(2).^2)./beta(2)^2;
score = trimmean(W.*abs(residuals).^2,20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [yhat,residuals] = fit(beta,x,y);
yhat = beta(1) + beta(2).*x;
residuals = (y - yhat);
