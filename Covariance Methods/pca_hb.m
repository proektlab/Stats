function [pc,z,latent]=pca(X)
%
% Calculates principal components
% Input: 
% X - data (samples x channels)
% Outputs:
% pc - principal component vectors
% z -  projections of data onto principal components
% latent - eigenvalues of the covariance matrix

C=covmat(X);
[pc,latent]=eig(C);
latent=diag(latent);
z=X*pc;
pc=pc(:,end:-1:1);
latent=latent(end:-1:1);
z=z(:,end:-1:1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C=covmat(X)
%
% Computes covariance matrix of high dimensional data after detrending
% Input:
% X - data (samples x channels)
% Output:
% C - covariance matrix

[N,Ch]=size(X);
if Ch==1; error('Need multidimensional data'); end
if N < Ch; error('Need more samples than channels'); end;
av=mean(X);
X=X-av(ones(N,1),:); % detrend
C=X.'*X;
% for i=1:Ch
%     for j=1:Ch
%         C(i,j)=X(:,i).'*X(:,j);
%     end;
% end;
C=C/(N-1);