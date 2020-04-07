function [m_hat,k_hat,m_ci,k_ci] = vonMisesFit(theta,alpha);
% function [m_hat,k_hat,m_ci,k_ci] = vonMisesFit(theta,alpha);


if nargin<2 alpha = []; end;
if isempty(alpha) alpha = 0.05; end;

n = size(theta,1);
C_bar = mean(cos(theta),1);
S_bar = mean(sin(theta),1);
theta_bar = atan2(S_bar,C_bar);
R2 = C_bar.^2+S_bar.^2;
R = sqrt(R2);
m_hat = theta_bar;
k_hat  = (1.28 - 0.53.*R2).*tan(pi*R./2);

m_ci = [theta_bar-acos((1-chi2cdf(1,1-alpha))./(2.*k_hat.*R)) theta_bar+acos((1-chi2cdf(1,1-alpha))./(2.*k_hat.*R))];

a = (n-R)./max(chi2cdf(1-alpha/2,n-1),eps);
b = (n-R)./max(chi2cdf(alpha/2,n-1),eps);
k_ci = [(1+(1+3*a).^0.5)./(4*a) (1+(1+3*b).^0.5)./(4*b)];
