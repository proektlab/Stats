function [output] = cca(X,Y);
% function [output] = cca(X,Y);
%
% This function performs the Canonical Correlation Analysis due to Harold
% Hotelling.  Each row of the data matrices X and Y is assumed to be a set
% of paired observations, so that size(X,1) = size(Y,1).
%
% Canonical Correlation Analysis searches for the set of orthogonal linear 
% transformations (a_j,b_j) that maximize the correlations, rho_j, between 
% the transformed variables. Thus, rho(j) = corr(X*a(:,j),Y*b(:,j)) is the
% maximum correlation possible for transformations (a(:,j), b(:,j)) that 
% is orthogonal to the previous j-1 transformations.
%
% This function, CCA, operates on the raw data.  If you have already
% computed the covariance matrices (or if you want to provide regularized
% covariance matrices, as in PDA), use the function
% ccaFromCovarianceMatrices, which this function calls.
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


[Nsamplesx,Mvariables] = size(X);
[Nsamples,Jvariables] = size(Y);

if Nsamples~=Nsamplesx,
    error('Canonical Correlation Analysis requires matched pairs of observations in the rows of X and Y.');
elseif Nsamples == 1
    error('Canonical Correlation Analysis requires more than one set of observations.');
end


grand_meanX = mean(X,1);
grand_meanY = mean(Y,1);

H = X - repmat(grand_meanX,[Nsamples 1]); % remove the mean response
I = Y - repmat(grand_meanY,[Nsamples 1]);

sig_11 = (1/Nsamples)*H.'*conj(H); 
sig_22 = (1/Nsamples)*I.'*conj(I); 
sig_12 = (1/Nsamples)*H.'*conj(I); % covariation of predictors and class variables
sig_21 = (1/Nsamples)*I.'*conj(H);

sig_11 = 0.5*(sig_11+sig_11'); % minimize roundoff assymmetry
sig_22 = 0.5*(sig_22+sig_22'); 
sig_12 = 0.5*(sig_12+sig_21');

output = ccaFromCovarianceMatrices(sig_11,sig_22,sig_12,Nsamples);
% Note: it is possible to do the entire CCA without explicitly calculating
% the correlation matrices via a set of svds.  I haven't compared the
% numerical stability of the approaches, though.