function [dz,vdz,Adz]=jackknifeTwoGroupCoherencies(coh1,coh2,p)
% function [dz,vdz,Adz]=jackknifeTwoGroupCoherencies(coh1,coh2,p)
%
% Inputs:
% J1c1   tapered fourier transform of dataset 1 in condition 1
% J2c1   tapered fourier transform of dataset 1 in condition 1
% J1c2   tapered fourier transform of dataset 1 in condition 2
% J2c2   tapered fourier transform of dataset 1 in condition 2
% p      p value for test (default: 0.05)
%
%
% Dimensions: J1c1,J2c2: frequencies x number of samples in condition 1
%              J1c2,J2c2: frequencies x number of samples in condition 2
%              number of samples = number of trials x number of tapers
% Outputs:
% dz    test statistic (will be distributed as N(0,1) under H0
% vdz   Arvesen estimate of the variance of dz
% Adz   1/0 for accept/reject null hypothesis of equal population
%       coherences based dz ~ N(0,1)
% 
% Note: all outputs are functions of frequency
%
% References: Arvesen, Jackkknifing U-statistics, Annals of Mathematical
% Statisitics, vol 40, no. 6, pg 2076-2100 (1969)

m1=size(coh1,1); % number of samples, condition 1
m2=size(coh2,1); % number of samples, condition 2
dof1=2*m1; % number of degrees of freedom in the first condition estimates
dof2=2*m2; % number of degrees of freedom in the second condition estimates

if nargin < 3; p=0.05; end; % set the default p value

Cm1=abs(mean(coh1,1)); % mean coherence, condition 1
Cm2=abs(mean(coh2,1)); % mean coherence, condition 2

%
% Compute the statistic dz, and the probability of observing the value dz
% given an N(0,1) distribution i.e. under the null hypothesis
%
z1=atanh(Cm1)-1/(dof1-2); % Bias-corrected Fisher z, condition 1
z2=atanh(Cm2)-1/(dof2-2); % Bias-corrected Fisher z, condition 2
dz=(z1-z2)/sqrt(1/(dof1-2)+1/(dof2-2)); % z statistic
%
% The remaining portion of the program computes Jackknife estimates of the mean (mdz) and variance (vdz) of dz
% 
samples1=[1:m1];
samples2=[1:m2];

coh1dropj = zeros(size(coh1));
z1i = zeros(size(coh1));
dz1i = z1i;
coh2dropj = zeros(size(coh2));
z2i = zeros(size(coh2));
dz2i = z2i;
%
% Leave one out of one sample
%
for i=1:m1;
    ikeep=setdiff(samples1,i); % all samples except i
    coh1dropj(i,:)=squeeze(abs(mean(coh1(ikeep,:),1))); % 1 drop mean spectrum, data 2, condition 1
    z1i(i,:,:)=atanh(coh1dropj(i,:,:))-1/(dof1-4); % 1 drop, bias-corrected Fisher z, condition 1
    dz1i(i,:,:)=(z1i(i,:,:)-z2)/sqrt(1/(dof1-4)+1/(dof2-2)); % 1 drop, z statistic, condition 1
    ps1(i,:,:)=m1*dz-(m1-1)*dz1i(i,:,:);
end; 
ps1m=mean(ps1,1);

for i=1:m2;
    ikeep=setdiff(samples2,i); % all samples except i
    coh2dropj(i,:)=squeeze(mean(coh2(ikeep,:),1)); % 1 drop mean spectrum, data 2, condition 1
    z2i(i,:,:)=atanh(coh2dropj(i,:,:))-1/(dof2-4); % 1 drop, bias-corrected Fisher z, condition 1
    dz2i(i,:,:)=(z1-z2i(i,:,:))/sqrt(1/(dof1-2)+1/(dof2-4)); % 1 drop, z statistic, condition 1
    ps2(i,:,:)=m2*dz-(m2-1)*dz2i(i,:,:);
end; 
ps2m=mean(ps2,1);

vdz=sum((ps1-ps1m(ones(m1,1),:,:)).*(ps1-ps1m(ones(m1,1),:,:)),1)/(m1*(m1-1))+sum((ps2-ps2m(ones(m2,1),:,:)).*(ps2-ps2m(ones(m2,1),:,:)),1)/(m2*(m2-1));

%
% Test whether H0 is accepted at the specified p value
%
Adz=zeros(size(dz));
x=norminv([p/2 1-p/2],0,1);
indx=find(dz>=x(1) & dz<=x(2)); 
Adz(indx)=1;
