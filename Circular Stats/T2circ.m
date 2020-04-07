function [p,radii,zhat,t2c] = T2circ(z,zeta,alpha);
% function [p,radii,zhat,t2c] = T2circ(z,zeta,alpha);
%
% T2circ calculates the statistic T^2_circ from Victor and Mast (1991).
% T^2_circ is designed to assess the variability of mulitple estimates of 
% Fourier components to assign significance to the estimates.  The
% fundamental assumption of T^2_circ is that the real and imaginary
% components of the Fourier components are independent and
% Gaussian-distributed with equal variances.  See the original reference
% for more detail:
%
%  Victor, J.D., and Mast, J. (1991) A new statistic for steady-state evoked 
%     potentials. Electroenceph. Clin. Neurophysiol. 78, 378-388. 
%
% INPUTS:
%     z, a vector or matrix of complex numbers (Fourier components). Each
%        row is assumed to be a separate estimate of the component.  If
%        multiple columns are present, t2c is calculated for each column
%        separately.
%     zeta, a scalar or vector, containing the expected value of z.
%        Defaults to 0.
%     alpha, a scalar, to determine the significance level for the
%        confidence limits. Defaults to 0.05.
%
% OUTPUTS:
%     p, a vector of probabilities that the population of estimates z could
%        come from a process with expected value zeta.  If p<alpha, then
%        the null hypothesis that E<z> == zeta can be rejected for that
%        probability level - that is, the closer to 0, the more significant.
%     radii, the radius of the alpha confidence interval around the point
%        zhat.  The 100*(1-alpha)% confidence interval of E<z(:,i)> is a 
%        circle centered on zhat(i) with a radius of radii(i).
%     zhat, the plug-in estimate of E<z>, namely, mean(z,1).
%     t2c, the T^2_circ statistic for each column of z; theoretically
%        distributed as F_2,_2*(# of estimates - 1).
%
% Function originally written by Andrew Hudson, 3/19/2005.


if nargin<2 zeta = []; end;
if nargin<3 alpha = []; end;

if isempty(zeta) zeta = 0; end;
if isempty(alpha) alpha = 0.05; end;

[M,N] = size(z);
if M == 1
    warning('T2circ: you passed a row vector in; will assume that these are estimates of a single variable');
    z = squeeze(z);
    [M,N] = size(z);
end

zeta = zeta(:)';
if ~(length(zeta==N)|length(zeta==1))
    error('T2circ: Zeta must either contain a scalar, or a vector with one entry for each column of z.');
end

if M<3 error('T2circ: Estimation of T^2_circ requires at least 3 estimates'); end;

if alpha>0.5, alpha = 1 - alpha; end;

zhat = mean(z,1);
residuals = z - repmat(zhat,[M 1]);
sumSqResid = sum(abs(residuals).^2,1);

t2c = (M-1)*abs(zhat - zeta).^2./sumSqResid;
p = 1 - fcdf(M*t2c,2,2*(M-1));
if nargout>1
    radii = sqrt(finv(1-alpha,2,2*(M -1))*sumSqResid/(M*(M-1)));
end
return;