function [jksd,jkHat] = powerSpectrumStandardErrorFromJackPop(jk);

N = size(jk,1);
jkHat = mean(10*log10(jk),1);
jksd = sqrt(((N-1)/N).*sum((10*log10(jk) - repmat(jkHat,[N 1 1])).^2,1));
