function [latency,latencySE] = angularRegress(data,cnum);

crf = extractCRF(data.cluster{1}.lineComponents,'Chat');
jcrf = extractCRF(data.cluster{1}.jkLineComponents,'Chat');

useHarmonics = [1 2 3 4 5];%1:5;
harmonics = (1:23).*96/46*2*pi;

for c = 1:6,
    if c>4, useHarmonics = [1 3 4 5]; end;
    [latency{1,1}(c),Dev(c),Rsq(c)] = computeLatency(harmonics(useHarmonics),angle(crf.allN11(23*(c-1)+useHarmonics))');
    include = find(jcrf.allC11==c);
    latencySE{1,1}(c) = jackknifeLatency(harmonics(useHarmonics),angle(jcrf.allN11(include,useHarmonics)));
end
for c = 1:6,
    if c>4, useHarmonics = [1 3 4 5]; end;
    [latency{2,1}(c),Dev(c),Rsq(c)] = computeLatency(harmonics(useHarmonics),angle(crf.allN21(23*(c-1)+useHarmonics))');
    include = find(jcrf.allC21==c);
    latencySE{2,1}(c) = jackknifeLatency(harmonics(useHarmonics),angle(jcrf.allN21(include,useHarmonics)));
end
for c = 1:6,
    if c>4, useHarmonics = [1 3 4 5]; end;
    [latency{1,2}(c),Dev(c),Rsq(c)] = computeLatency(harmonics(useHarmonics),angle(crf.allN12(23*(c-1)+useHarmonics))');
    include = find(jcrf.allC12==c);
    latencySE{1,2}(c) = jackknifeLatency(harmonics(useHarmonics),angle(jcrf.allN12(include,useHarmonics)));
end
for c = 1:6,
    if c>4, useHarmonics = [1 3 4 5]; end;
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
 w = 1-dev;
 latencySE = sqrt((N-1)/N*sum(w.*(late - late*w'./sum(w)).^2)./sum(w));
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [late,dev,rsq] = computeLatency(X,Y);
 %X = X - X(1);
 X = X - mean(X);
 if mod(length(X),2)==1
     Y = Y - Y(ceil(length(X)/2));
 else
     Y = Y - mean(Y(length(X)/2:length(X)/2+1));
 end
 keyboard;
 seed = [0 late2param(mean(mod(diff(Y),2*pi))./diff(X(1:2)))];
 for j=1:25,  % multiple repeats to avoid local minima
    params{j} = nlinfit(X,Y,@ar,seed+2*(rand(1,2)-0.5));
    err(j) = deviance(ar(params{j},X),Y);
 end;
 late = param2late(params{min(find(err == nanmin(err)))}(2));
 dev = err(min(find(err==nanmin(err))));
 predict = ar(params{min(find(err==min(err')))},X);
 rsq = predict*Y'/sqrt(predict*predict')/sqrt(Y*Y');
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function err = deviance(X,Y);
 err = sqrt(mean(sin(abs(Y-X)/2).^2));
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theta = ar(params,X);
 params = real(params);
 %m = pi*(tanh(params(1)));
 m = params(1);
 b = param2late(params(2));
 theta = angle(exp(i*(X*b+m)));
 %theta = (exp(i*(X*b+m)));
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function late = param2late(param);
 param = real(param);
 late = (46/96*(1+tanh(param/10))/2);
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function param = late2param(late);
 late = real(late);
 param = 10*atanh(96*2/46*late-1);
return;