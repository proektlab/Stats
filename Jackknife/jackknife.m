function [muHat,sigmaHat]=jackknife(data);
%function [muHat,sigmaHat]=jackknife(data);
%
% This function computes the mean and standard error of a dataset.
% Each row of the data is assumed to be an independent trial.  The data may
% be a vector or multidimensional.
%
% Created by Andrew Hudson, 11/21/2003.

nd = ndims(data);
dims=size(data);

N = dims(1);
% Create jackknife population
if length(dims)>2
   jk = reshape((1-eye(N))*reshape(data, [N prod(dims(2:end))])/(N-1),dims);
else
   jk = (1-eye(N))*data./(N-1);
end
% Make relevant calculations
muHat=mean(jk,1);
sigmaHat=sqrt(((N-1)/N).*sum((jk - repmat(muHat,[N ones(1, nd-1)])).^2,1));
biasHat = (N-1)*(mean(data,1) - muHat);