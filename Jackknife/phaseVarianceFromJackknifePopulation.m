function [jkCircStandardErr,jkCircVariance] = phaseVarianceFromJackknifePopulation(jkpop);
% function [jkCircStandardErr,jkCircVariance] = phaseVarianceFromJackknifePopulation(jkpop);
%
% INPUT:  jkpop, a [# of trials, # of variables] matrix of jackknife
%               estimates of the population mean
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
% This function should be used when a population of jackknife estimates is
% already available.  If instead you have a set of raw data to jackknife,
% use the function jaccknifePhaseVariance.
%
% Written by Andrew Hudson (c), 7/8/2006.
if isreal(jkpop), jkpop = exp(i*jkpop); end;


N = size(jkpop,1);

jkCircVariance = (N-1).*circularVariance(jkpop); % This is the jackknife
   % version of the circular variance based upon Marida and Jupp (2000), 
   % Directional Statistics.  Note that Thomson and Chave (1991) 
   % "Jackknifed Error Estimates for Spectra, Coherences, and Transfer
   % Functions" in Advances in Spectrum Analysis and Array Processing, 
   % Vol 1, S Haykin, Ed. define the circular variance as 2*circVariance.
   
jkCircStandardErr = zeros(size(jkCircVariance));
jkCircStandardErr(find(jkCircVariance<1)) = sqrt(-2*log(1-jkCircVariance(find(jkCircVariance<1)))); 
  % approx = sqrt(2*circVariance) for small circVariance.
  % Circular standard deviation is useful for plotting phase errors
  % on the line, while circular variance is more useful on the circle.

% Now correct for the spread added by the jackknife estimation, which can result in complex numbers  
% the circular standard deviation can dake values in [0,+Inf).
% circVariance takes values in [0,1].
jkCircStandardErr(find(jkCircVariance>=1)) = Inf;
jkCircVariance = min(jkCircVariance,1);