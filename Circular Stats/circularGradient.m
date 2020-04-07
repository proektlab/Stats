function grad = circularGradient(theta);
% function grad = circularGradient(theta);
%
% Uses the circular distance to determine an estimate of the gradient of
% theta, along the first non-singleton dimension of theta.  The estimate
% grad is determined using the forward and reverse single difference at the
% ends of the vector, and taking the circular average of the single
% differences before and after interior points.
% 
% Written by Andrew Hudson, (c) 7/11/2006.

nd = ndims(theta);
dims = size(theta);
along = min(find(dims>1));

theta = shiftdim(theta,along-2);

dt = circularDistance(theta(1,1:end-1,:),theta(1,2:end,:));
dtm = circularMean([dt(1,1:end-1);dt(1,2:end)]);
grad = shiftdim([dt(1);dtm';dt(end)],2);