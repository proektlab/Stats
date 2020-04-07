function [latency,latencySE] = angularRegress(data,cnum,squarewave);

crf = extractCRF(data.cluster{1}.lineComponents,'Chat');
jcrf = extractCRF(data.cluster{1}.jkLineComponents,'Chat');

useHarmonics = [1 2 3 4 5];%1:5;
harmonics = (1:23).*96/46*2*pi;

if nargin<3
    squarewave.cs = ones(1,23);
end
squarewave.cs = squarewave.cs(useHarmonics);

for c = 1:6,
    [latency{1,1}(c),Dev(c),Rsq(c)] = computeLatency(harmonics(useHarmonics),angle(crf.allN11(23*(c-1)+useHarmonics)'.*conj(squarewave.cs)./abs(squarewave.cs).^2));
    include = find(jcrf.allC11==c);
    latencySE{1,1}(c) = jackknifeLatency(harmonics(useHarmonics),angle(jcrf.allN11(include,useHarmonics).*repmat(conj(squarewave.cs)./abs(squarewave.cs).^2,[length(include) 1])));
end
keyboard;
for c = 1:6,
    [latency{2,1}(c),Dev(c),Rsq(c)] = computeLatency(harmonics(useHarmonics),angle(crf.allN21(23*(c-1)+useHarmonics))');
    include = find(jcrf.allC21==c);
    latencySE{2,1}(c) = jackknifeLatency(harmonics(useHarmonics),angle(jcrf.allN21(include,useHarmonics)));
end
keyboard;
for c = 1:6,
    [latency{1,2}(c),Dev(c),Rsq(c)] = computeLatency(harmonics(useHarmonics),angle(crf.allN12(23*(c-1)+useHarmonics))');
    include = find(jcrf.allC12==c);
    latencySE{1,2}(c) = jackknifeLatency(harmonics(useHarmonics),angle(jcrf.allN12(include,useHarmonics)));
end
for c = 1:6,
    [latency{2,2}(c),Dev(c),Rsq(c)] = computeLatency(harmonics(useHarmonics),angle(crf.allN22(23*(c-1)+useHarmonics))');
    include = find(jcrf.allC22==c);
    latencySE{2,2}(c) = jackknifeLatency(harmonics(useHarmonics),angle(jcrf.allN22(include,useHarmonics)));
end
keyboard;

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function latencySE = jackknifeLatency(X,crf);
 N = size(crf,1);
 for i=1:N
    Y = crf(i,:);
    [late(i),dev(i),rsq(i)] = computeLatency(X,Y);
 end
 latencySE = sqrt((N-1)/N*sum(late - mean(late)).^2);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [late,dev,rsq] = computeLatency(X,Y);
 

 N = length(X);
 if Y(1) > 0, Y(1) = Y(1) - 2*pi; end;
 Y = Y(1)+rewrap(angle(exp(i.*Y).*exp(-i*Y(1))));
 % Y = rewrap(Y);
 % X = [-X(end:-1:1)  X];
 % Y = [-Y(end:-1:1)  Y];
 X = [zeros(size(X,1),1)  X];
 Y = [zeros(size(Y,1),1)  Y];
 
 %late = regress(Y',X');
 %late = regress(Y',[ones(size(X')) X']);
 %predict = X'*late;
 %late = -late;
 late = robustfit(X,Y);
 predict = [ones(size(X')) X']*late;
 late = -late(2);
 dev = deviance(predict,Y);
 rsq = circcorr(predict,Y);
 
%  seed = [late2param(mean(mod(diff(Y),2*pi))./diff(X(1:2)))];
%  for j=1:25,  % multiple repeats to avoid local minimab
%     params{j} = fminsearch(@score,seed+2*(rand(1,1)-0.5),[],X,Y);
%     % params{j} = nlinfit(X,Y,@ar,seed+2*(rand(1,1)-0.5));
%     err(j) = deviance(ar(params{j},X),Y);
%  end;
%  late = param2late(params{min(find(err == nanmin(err)))});
%  dev = err(min(find(err==nanmin(err))));
%  predict = ar(params{min(find(err==min(err')))},X(N+1:end));

%  dev = deviance(predict(N+1:end),Y(N+1:end));
%  rsq = circcorr(predict(N+1:end),Y(N+1:end));
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Yo = rewrap(Y);
jump = -0.52;
Yo = Y;
DYo = diff([zeros(size(Y,1),1) Y],1,2);
while(any(DYo>jump)),
    Yo = Yo - 2*pi*cumsum(DYo>jump,2);
    DYo = diff([zeros(size(Y,1),1) Yo],1,2);
end;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = score(param,X,Y);
 W = [1:length(X)/2, length(X)/2:-1:1];
 W = sqrt(W)./max(sqrt(W));
 err = deviance(Y,ar(param,X),W);
  %err = 2 - circcorr(Y,ar(param,X));
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhosq = circcorr(X,Y);
 % From Jupp and Mardia, 1980. (Biometrika 67:1, 163-173).
 Cx = cos(X);
 Sx = sin(X);
 Cy = cos(Y);
 Sy = sin(Y);
 
 Rcc = correlation(Cx,Cy);
 Rcs = correlation(Cx,Sy);
 Rsc = correlation(Sx,Cy);
 Rss = correlation(Sx,Sy);
 
 R1 = correlation(Cx,Sx);
 R2 = correlation(Cy,Sy);
 
 
 rhosq = (Rcc.^2+Rcs.^2+Rsc.^2+Rss.^2 + 2*(Rcc*Rss - Rcs*Rsc)*R1*R2 ...
     - 2*(Rcc*Rcs + Rsc*Rss)*R2 - 2*(Rcc*Rsc + Rcs*Rss)*R1) ...
     ./ ((1 - R1.^2)*(1 - R2.^2));
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = correlation(X,Y);
 X = X(:); Y = Y(:);
 r = corrcoef(X,Y);
 r = r(1,2);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = circDistance(X,Y);
 X = X(:); Y = Y(:);
 d = 1 - cos(X-Y);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = deviance(X,Y,w);
 if nargin<3,w = []; end;
 X = X(:); Y = Y(:);
 if isempty(w), w = ones(1,length(X)); end;
 err = w*circDistance(X,Y)/sum(w);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta = ar(params,X);
 params = real(params);
 b = param2late(params(1));
 theta = angle(exp(i*(-X*b)));
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function late = param2late(param);
 param = real(param);
 late = (46/96/2*(1+tanh(param/10))/2);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function param = late2param(late);
 late = real(late);
 param = 10*atanh(96*4/46*(late)-1);
return;