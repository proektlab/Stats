function [fs,cs,sp]=ftest(ts,nw,k,npad)

ts = detrend(ts(:)); % remove any linear trends in the data
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



xk = taperedSpectralEstimate(ts',v,npad,0.001)/sqrt(2*0.001);
% xk1 is the set of tapered spectral estimates
% arranged as [channel,frequency,taper]
Hk0 = sum(v(:,1:2:k)); % what is up with taking every other taper?
Sk = abs(xk.^2);

Sp = squeeze(sum(Sk,3));
for ch = 1:size(ts,2)
    Cs(ch,:) = (squeeze(xk(ch,:,1:2:k))*Hk0'/sum(Hk0.^2))';
    Num = abs(Cs(ch,:).^2);
    Fs(ch,:) = (k-1)*Num./(Sp(ch,:)/sum(Hk0.^2)-Num);
end
keyboard;
