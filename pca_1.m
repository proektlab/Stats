function [T,pvar,W,L] = pca_1(Y)
%PCA Principal Component Analysis (PCA)
%   T = PCA(Y) returns the principal components T (i.e., temporal
%   eigenvectors) of the matrix Y.
%
%   [T,PVAR,W,L] = PCA(Y) returns optional output according to the
%   following equations: Y = C'WT and T = L*Y.  Where C are the spatial
%   eigenvectors, W is a diagonal matrix of the eigenvalues (the
%   variance in each component), and L = inv(W)*C is the linear
%   tranformation used to obtain the temporal eigenvectors T from
%   the matrix Y. PVAR is the percent variance explained by each
%   component.
%
%   According to PCA, C and T are orthonormal (i.e., independent or
%   instantaneously uncorrelated and scaled to unity).
%
%   See also SVD.
%
%   Michael A. Repucci
%   Last Modified: 5/25/2000 10:57 PM

[c,t] = size(Y);
%if c>t
%   Y = Y';
%   tmp = c;
%   c = t;
%   t = tmp;
%end

% The singular value decomposition (SVD) is the backbone of PCA
[U,SS,V] = svd(Y*Y');
W = sqrt(SS);
w = diag(W);
L = diag(1./w)*U';
T = L*Y;

% The percent variance (PVAR) explained by each principle component is
% a fraction of the total variance
pvar = 100*(w.^2)/sum(w.^2);
