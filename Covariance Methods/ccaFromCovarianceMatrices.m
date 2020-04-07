function [output,isig_11,isig_22] = ccaFromCovarianceMatrices(sig_11,sig_22,sig_12,Nsamples);
% function output = ccaFromCovarianceMatrices(sig_xx,sig_yy,sig_xy,Nsamples); 
%
% This function computes the canonical correlations between two sets of
% variables from their covariance matrices (sig_xx, sig_yy), their cross-
% covariance matrix, sig_xy, and the number of repeated observations,
% Nsamples.
%
% The results of the canonical correlation analysis are returned in the
% structure output, with the following fields:
%
%    rho, the canonical correlations. The fraction of explained
%         cross-covariance is rho.^2
%    a, a matrix containing the transforms to convert the first set of
%         variables into canonical variates.  Each set of transforms is
%         arranged in a column of a, so the ith canonical transform of 
%         variable x is x*a(:,i).  
%    b, a matrix containing the transforms to convert the second set of
%         variables into canonical variates.  Each set of transforms is
%         arranged in a column of b, so the ith canonical transform of 
%         variable y is y*b(:,i).  
%    r2a,r2b; vectors indicating the fraction of explained covariance 
%         captured by the canonical variates (x*a(:,i), y*b(:,i)).  Hence,
%         r2a(k) indicates the fraction of sig_xx captured by the kth
%         canonical variate.
%    The remaining field is only calcluated if Nsamples is provided:
%    P, The probability that the first k canonical correlations are
%         nonzero.  Hence, values close to one suggest the first k
%         canonical correlations are nonzero.
%
% Written by Andrew Hudson (c) 7/17/2006.
if nargin<4
    Nsamples  = [];
end

Mvariables = size(sig_11,1);
Jvariables = size(sig_22,1);

retain = min(Mvariables,Jvariables);

% find the inverse square roots of the correlation matrices
[eig11,lambda11] = eig(sig_11); 
iss11 = eig11*inv(sqrt(lambda11))*inv(eig11);
iss11 = 0.5*(iss11+iss11'); 

[eig22,lambda22] = eig(sig_22); 
iss22 = eig22*inv(sqrt(lambda22))*inv(eig22);
iss22 = 0.5*(iss22+iss22'); 

if nargout>1
    isig_11 = eig11*inv(lambda11)*inv(eig11);
    isig_22 = eig22*inv(lambda22)*inv(eig22);
end

% Solve the CCA problem
[t_1,Da,b_1] = svd(iss11*sig_12*iss22);

a = iss11*t_1;
b = iss22*b_1;

ia = inv(a.');
ib = inv(b.');

a = a(:,1:retain);
b = b(:,1:retain);

rho = diag(Da);
rho = rho(1:retain);

% Diagnostics
% Proportion of standardized sample variance explained by each canonical variate
r2a = diag(ia'*ia)./trace(ia'*ia); 
r2b = diag(ib'*ib)./trace(ib'*ib);

output.rho = rho;
output.a = a;
output.b = b;
output.r2a = r2a;
output.r2b = r2b;

if ~isempty(Nsamples)
    % Likelihood ratio test for non-zero canonical correlations
    for k = 1:retain % test to see whether the first k cc's are non-zero
        chi2k(k) = -(Nsamples - 1 - 0.5*(Mvariables + Jvariables +1))*sum(log(1 - rho(k:retain).^2));
        nuChi2k(k) = (Mvariables-k+1)*(Jvariables-k+1);
    end;
    output.P = chi2cdf(chi2k,nuChi2k);
else
    output.P = [];
end
