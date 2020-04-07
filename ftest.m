function [fs,cs,sp]=ftest(ts,nw,k,npad)

% Compute F-statistic for sine wave in locally white noise
% model
ts=ts(:);
nt=length(ts);

v=dpss(nt,nw,k);
ind=[1:2:k];
wk0=sum(v(:,ind));
ts1=repmat(ts,1,k);
tsft=fft(ts1.*v,npad);
sp=sum(abs(tsft.^2)')';
cs=tsft(:,ind)*wk0'/sum(wk0.^2);
num=abs(cs.^2);
fs=num./(sp/sum(wk0.^2)-num)*(k-1);
