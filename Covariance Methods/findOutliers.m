function [outliers,bfe] = findOutliers(data,tolerance,pccrit,criterion);
% function [outliers,bfe] = findOutliers(data,tolerance,pccrit,criterion);
%
% Ok, the idea behind this approach is to keep dropping outliers until the 
% principal components of the covariance matrix stop shifting their orientation.
% Then these components are used to compute tsquared for all of the points, and
% those that are above a certain threshold will be declared outliers.

if nargin<2
   tolerance = pi/180; % start with 1 degree
end
if nargin<3
   pccrit = 0.99; % for the iterating stage
end
if nargin<4
   criterion = 0.99995; % for the exclusion stage
end
outlierMultiple = 1.5; % typically 1.5 or 3
% seed with initial values
[nSamples,nVars]=size(data);
[pc,~,lambda,tsq] = princomp(data);
p=fcdf(tsq,nVars,nSamples-nVars);
projection = data*pc;
iqrp = iqr(projection);
prctiles = prctile(projection,[25 75]);
limits = [prctiles(1,:)-outlierMultiple*iqrp;prctiles(2,:)+outlierMultiple*iqrp];

% start by flagging everything that's a long way out
outliers=(find(or(or(projection(:,1)<limits(ones(nSamples,1),1),projection(:,1)>limits(2*ones(nSamples,1),1)),...
   or(projection(:,2)<limits(ones(nSamples,1),2),projection(:,2)>limits(2*ones(nSamples,1),2)))));

samples = 1:nSamples;
%outliers = samples(find(p>pccrit));
opc = pc; % old PCs
olambda = lambda;

done = 0;
while ~done
   % drop outliers
   subset = setdiff(samples,outliers);
   N = length(subset);
   [npc,~,nlambda,ntsq] = princomp(data(subset,:));
   np=fcdf(ntsq,nVars,N-nVars);
   noutliers = subset(np>pccrit);
   lopc = sqrt(opc(:,1)'*opc(:,1));
   lnpc = sqrt(npc(:,1)'*npc(:,1));
   pcangle = acos(npc(:,1)'*opc(:,1)./(lopc*lnpc)); % the angle between the principle components
   disp(['Angle is ' num2str(pcangle)]);
   if pcangle > tolerance % then save the new values and loop again
      outliers = union(outliers,noutliers);
      opc = npc;
      olambda = nlambda;
   else
      done = 1;
   end
end
subset = setdiff(samples,outliers);

% now we have found a stable set of eigenvectors and we can use these to find outliers
avgData = mean(data(subset,:));
centeredData = data - avgData(ones(nSamples,1),:);
score = centeredData*opc; % project onto the principal components
tmp =  sqrt(diag(1./olambda))*score';
alltsq = sum(tmp.*tmp)';
allp=fcdf(alltsq,nVars,nSamples-nVars);
outliers = samples(allp>criterion);
if nargout > 1
    [~,bfe]=t2confidenceRegion(data(subset,:),criterion);
end;