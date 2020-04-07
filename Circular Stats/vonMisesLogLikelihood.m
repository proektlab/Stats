function l = vonMisesLogLikelihood(theta,M,K);
% function l = vonMisesLogLikelihood(theta,M,K);
%

if size(theta,2) == length(M)
    M = repmat(M(:)',[size(theta,1) 1]);
end
if size(theta,2) == length(K)
    K = repmat(K(:)',[size(theta,1) 1]);
end

l = K.*cos(theta-M) - log(2*pi.*besseli(0,K));