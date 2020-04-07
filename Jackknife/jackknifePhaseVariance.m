function [jkCircStandardErr,jkCircVariance] = jaccknifePhaseVariance(data);
% function [jkCircStandardErr,jkCircVariance] = jaccknifePhaseVariance(data);
%
% INPUT:  data, a [# of trials, # of variables] matrix of complex-valued
%               data
% OUTPUT: jkCircStandardErr, a (# of variables) element vector containing
%               an estimate of spread that ranges from [0,+Inf), for
%               plotting on the real line.
%         jkCircVariance, a (# of variables) element vector containing an
%               estimate of spread that ranges from [0,1], for plotting on
%               the circle.  A value of 1 indicates that the data are
%               uniformly spread about the circle.
%
% This function provides an estimate of the phase spread in the population.
% Since this function was designed with spectral estimates in mind, the
% data are assumed to live in the complex plane.
%
% This function should be used when you have a set of raw data to jackknife.  
% If instead you have a population of jackknife estimates already, use the 
% function phaseVarianceFromJackknifePopulation.
%
% Written by Andrew Hudson (c), 7/8/2006.


dims=size(data);

N = dims(1);
% Create jackknife population
if length(dims)>2
   jk = reshape((1-eye(N))*reshape(data, [N prod(dims(2:end))])/(N-1),dims);
else
   jk = (1-eye(N))*data./(N-1);
end

[jkCircStandardErr,jkCircVariance] = phaseVarianceFromJackknifePopulation(jk);