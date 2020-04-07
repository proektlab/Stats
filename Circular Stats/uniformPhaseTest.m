function [significantPoints,significantValues,P]=uniformPhaseTest(statistic,criterion);
%function [significantPoints,significantValues,P]=uniformPhaseTest(statistic,criterion);
%
% This function will compare the values of a complex statistic against
%  a uniform phase distribution.
%
% INPUTS:
%  statistic: a complex-valued matrix.  The first dimension 
%             of this matrix is assumed to be multiple observations
%             thus a typical format for a spectrogram would be
%             [samples timepoints frequencypoints]
%  criterion: a probability criterion, e.g., 0.05, or 0.01, or a vector of criteria
%             defaults to 0.05
%
% Written by Andrew Hudson, 8/12/2003.  Last updated 6/1/2004.

dimensions = size(statistic);
num_dimensions = length(dimensions);
samples = dimensions(1);
if nargin<2
   criterion = 0.05;
end;

% first we make sure that we are looking at phase factors - i.e., that our
% statistic has a magnitude of 1 around the circle
nz = find(statistic);
statistic(nz) = statistic(nz)./abs(statistic(nz));

disp('Looking for significant deviations from a uniform phase distribution');
% Since the statistic is complex, taking the mean before going further estimates
% the position relative to the origin, which be the mean for a uniform phase distribution
Amean=squeeze(abs(mean(statistic,1))); % the norm of the mean statistic value
Pmean=squeeze(angle(mean(statistic,1))); % the phase of the mean statistic value

significantPoints=zeros([length(criterion) dimensions(2:end)]);

z = samples .*(Amean.^2); % Calculate the critical values of the Rayleigh distribution using the saddle-node approximation from Mardia & Jupp, p95
P = exp(-z).*(1+((2.*z) - z.^2)/(4.*samples) - (24*z - 132*z.^2+76*z.^3- 9*z.^4)/(288.*samples.^2));
for crit=1:length(criterion)
   if num_dimensions == 2
      what_signif = zeros(1,dimensions(2));
   else
      what_signif = zeros(dimensions(2:end));
   end
   if length(find(P<=criterion(crit)))
      what_signif(find(P<=criterion(crit)))=1;
   end
   significantPoints(crit,:)=what_signif(:)';
end

if length(dimensions)>2
    significantValues=significantPoints.*repmat(shiftdim(Amean,-1),[length(criterion) ones([1 length(dimensions(2:end))])]);
else
    significantValues=significantPoints.*Amean;
end