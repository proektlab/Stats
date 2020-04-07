function Sn = MADn(X,model,varargin);
% function Sn = MADn(X,model,varargin);
% 
% Calculates the median absolute deviation for a robust estimate of spread.
if nargin<2, model = 'norm'; end;

bn = fcorfac(X);
M = median(X,1);
d = modelScale(model,varargin{:});
Sn = bn*median(abs(X - repmat(M,[size(X,1) 1])),1)./d;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = modelScale(model,varargin);
if strcmp(lower(model),'norm')|strcmp(lower(model),'normal')|strcmp(lower(model),'gauss')
    med = norminv(0.5,varargin{:});
    d = fminsearch(@normMAD, norminv(.75,varargin{:}),[],med,varargin{:});
elseif strcmp(lower(model),'chi2')|strcmp(lower(model),'chisquare')
    med = chi2inv(0.5,varargin{:});
    d = fminsearch(@chi2MAD, chi2inv(.75,varargin{:}),[],med,varargin{:});
elseif strcmp(lower(model),'rayleigh')|strcmp(lower(model),'rayl')
    med = raylinv(0.5,varargin{:});
    d = fminsearch(@raylMAD, raylinv(.75,varargin{:}),[],med,varargin{:});
elseif strcmp(lower(model),'f')
    med = finv(0.5,varargin{:});
    d = fminsearch(@fMAD, finv(.75,varargin{:}),[],med,varargin{:});
else
    error('Unrecognized model');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = normMAD(x,med,varargin)
d = (normcdf(x+med,varargin{:}) - normcdf(x-med,varargin{:}) - 0.5).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = chi2MAD(x,med,varargin)
d = (chi2cdf(x+med,varargin{:}) - chi2cdf(x-med,varargin{:}) - 0.5).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = raylMAD(x,med,varargin)
d = (raylcdf(x+med,varargin{:}) - raylcdf(x-med,varargin{:}) - 0.5).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = fMAD(x,med,varargin)
d = (fcdf(x+med,varargin{:}) - fcdf(x-med,varargin{:}) - 0.5).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bn=fcorfac(Z)
n=size(Z,1);
switch n
case 2
   bn=1.196;
case 3
   bn=1.495;
case 4
   bn=1.363;
case 5
   bn=1.206;
case 6
   bn=1.200;
case 7
   bn=1.140;
case 8
   bn=1.129;
case 9
   bn=1.107;
end
if n>9
   bn=n/(n-0.8);
end