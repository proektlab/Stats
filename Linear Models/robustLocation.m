function location = robustLocation(X,modelname,varargin)
%function location = robustLocation(X,modelname,varargin)
%
% Returns a robust estimate of the central location of the data in X.
% possible models, specified by the string model name, include:
%       'norm','normal', or 'gauss': normal distribution
%       'chi2','chisquare': chi squared distribution
%       'rayleigh','rayl': rayleigh distribution
%       'f': F-distribution
% Additional arguments required to specify distributions (e.g., degrees
%    of freedom for F distribution) may be passed as vargin
% 
if nargin<2, model.name = []; end;
if isempty(modelname), model.name = 'norm'; 
else model.name = modelname; end;
model.d = [];
model.T1 = [];


nObs = size(X,1);
tol = .01;
maxIter = 500;

if nObs == 1,
    location = X;
elseif nObs == 2,
    location = mean(X);
elseif nObs ==3,
    location = median(X);
else
    T1 = median(X);
    model.T1 = T1;
    Sn = MADn(X,model,varargin{:});
    weightDist = psiPrimeDPhi(model,varargin{:});
    notConverged = 1;
    iter = 0;
    while(notConverged),
        iter = iter+1;
        T0 = T1;
        T1 = T0 + Sn.*mean(psiLog((X-repmat(T0,[nObs 1]))./repmat(Sn,[nObs 1])))./weightDist;
        notConverged = any((T1-T0)*(T1-T0)'>tol*T0*T0') && iter<maxIter;
    end
    if iter>=maxIter, warning('robustLocation did not converge; returning most recent results'); end;
    location = T0;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function psi = psiLog(X)
eX = exp(X);
psi = (eX-1)./(eX + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = psiPrimeLog(X)
y = 2.*exp(X)./(1+exp(X)).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Sn = MADn(X,model,varargin)
bn = fcorfac(X);
M = model.T1;
model = modelScale(model,varargin{:});
Sn = bn*median(abs(X - repmat(M,[size(X,1) 1])),1)./model.d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = modelScale(model,varargin)
if isempty(model.d)
    if strcmpi(model.name,'norm')||strcmpi(model.name,'normal')||strcmpi(model.name,'gauss')
        med = norminv(0.5,varargin{:});
        model.d = fminsearch(@normMAD, norminv(.75,varargin{:}),[],med,varargin{:});
    elseif strcmpi(model.name,'chi2')||strcmpi(model.name,'chisquare')
        med = chiinv(0.5,varargin{:});
        model.d = fminsearch(@chi2MAD, chi2inv(.75,varargin{:}),[],med,varargin{:});
    elseif strcmpi(model.name,'rayleigh')||strcmpi(model.name,'rayl')
        if varargin{:} == 1
            model.d = 2.351478951553029;
        else
            phi = raylinv(.75,varargin{:});
            med = raylinv(0.5,varargin{:});
            model.d = fminsearch(@raylMAD,phi,[],med,varargin{:});
        end
    elseif strcmpi(model.name,'f')
        if varargin{1} == 2 && varargin{2} == 2
            model.d = 1.236035156249995825561427;
        elseif varargin{1} == 2 && varargin{2} == 4
            model.d = 1.129687499999999289457264;
        elseif varargin{1} == 2 && varargin{2} == 6
            model.d = 1.110153570189508442922488;
        elseif varargin{1} == 2 && varargin{2} == 8
            model.d = 1.103733521768775815985464;
        elseif varargin{1} == 2 && varargin{2} == 10
            model.d = 1.10096060855679711565358;
        elseif varargin{1} == 2 && varargin{2} == 12
            model.d = 1.099587879047448257807673;
        else
            med = finv(0.5,varargin{:});
            phi = finv(0.75,varargin{:});
            model.d = fminsearch(@fMAD,phi,[],med,varargin{:});
        end
    else
        error('Unrecognized model.name');
    end;
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
function d = psiPrimeDPhi(model,varargin)
grid = 0.001:0.0025:0.999;
model.name = lower(model.name);
if strcmp(model.name,'norm')||strcmp(model.name,'normal')||strcmp(model.name,'gauss')
    grid = norminv(grid,varargin{:});
    d = normpdf(grid,varargin{:})*psiPrimeLog(grid)';
elseif strcmp(model.name,'chi2')||strcmp(model.name,'chisquare')
    grid = chi2inv(grid,varargin{:});
    d = chi2pdf(grid,varargin{:})*psiPrimeLog(grid)';
elseif strcmp(model.name,'rayleigh')||strcmp(model.name,'rayl')
    grid = raylinv(grid,varargin{:});
    d = raylpdf(grid,varargin{:})*psiPrimeLog(grid)';
elseif strcmp(model.name,'f')
    grid = finv(grid,varargin{:});
    d = fpdf(grid,varargin{:})*psiPrimeLog(grid)';
else
    error('Unrecognized model.name');
end;


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