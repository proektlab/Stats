function circVar = circularVariance(data);
% function circVar = circularVariance(data);
%
% INPUT:  data, a [# of trials, # of variables] matrix of complex data; to
%                use this function with data in the cartesian plane, first
%                create a data vector so that data = x+iy.
%
% OUTPUT: circVar, a (# of variables) element vector containing an
%               estimate of spread that ranges from [0,1], for plotting on
%               the circle.  A value of 1 indicates that the data are
%               uniformly spread about the circle.
%
% This function provides an estimate of the phase spread in the population.
% Since this function was designed with spectral estimates in mind, the
% data are assumed to live in the complex plane.
%
% This function should be used when a population of jackknife estimates is
% already available.  If instead you have a set of raw data to jackknife,
% use the function jaccknifePhaseVariance.
%
% Written by Andrew Hudson (c), 7/8/2006.
if all(isreal(data)), data = exp(i*data); end;
nz = find(data);
e = ones(size(data))*NaN;
e(nz) = data(nz)./abs(data(nz)); % convert to phase factors
e_resultant = mean(e,1);

circVar = (1 - abs(e_resultant));    % circVariance takes values in [0,1].
   % This is based upon Marida and Jupp (2000), Directional Statistics.  
   % Note that Thomson and Chave (1991) "Jackknifed Error Estimates for 
   % Spectra, Coherences, and Transfer Functions" in Advances in Spectrum 
   % Analysis and Array Processing, Vol 1, S Haykin, Ed. define the 
   % circular variance as 2*circVariance.