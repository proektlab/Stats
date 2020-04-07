function [ratem,rateci,timebins,rate] = rateEstimation(aligned,duration,pretime,binsize,nSpikes);
%function [ratem,rateci,timepts,rate] = rateEstimation(aligned,duration,pretime,binsize,nSpikes);
%
% INPUTS:
%     aligned, a cell array of spike times, one cell per trial
%     duration, a scalar, specifying the total amount of time in each entry of
%           the aligned structure in seconds
%     pretime, a scalar, indicating the amount of "negative time" included
%           in the aligned structure, in seconds
%     binsize, a scalar, specifying the time in seconds between the
%           interpolated rate values
%     nSpikes, a parameter that governs the bandwidth of the smoothing; the
%           nearest neighbor bandwidth will include this many spikes
%
% OUTPUTS:
%     ratem, the mean Poisson firing rate
%     rateci, the jackknife confidence intervals for the Poisson firing
%           rate
%     timepts, the time grid corresponding to each entry in ratem
%     rate, the jackknife population of rate estimates
if nargin<2, duration = []; end;
if nargin<3, pretime = []; end;
if nargin<4, binsize = []; end;
if nargin<5, nSpikes = []; end;

if isempty(duration), duration = max(cat(1,aligned{:})); end;
if isempty(pretime), pretime = min(cat(1,aligned{:})); end;
if isempty(binsize), binsize = 0.001; end;
if isempty(nSpikes), nSpikes = 100; end;

doLog = 0;
doSqrt = 1; % note: sqrt is the variance-stabilizing link for a poisson distribution
if iscell(aligned)
    Ntrials = length(aligned);
else
    Ntrials = 1;
    a{1} = aligned;
    aligned = a;
    clear a;
end

timebins = (binsize.*(1:round(duration./binsize))-pretime)';
if Ntrials>3
    wh=waitbar(0,'Local Rate Estimation');
    for j=1:Ntrials 
        waitbar(j/Ntrials,wh);
        spikes = sort(cat(1,aligned{[1:j-1 j+1:Ntrials]}));
        spikes = [-pretime - 1; spikes; duration-pretime+1];
        %[fit] = locfit(spikes,'alpha',[nSpikes/length(spikes) 0.05*0.75/2],'deg',2,'family','rate','mg',timebins,'xlim',[-pretime duration],'maxk',1000);
        %[fit] = locfit(spikes,'alpha',[nSpikes/length(spikes) 0.05*0.75/2],'deg',2,'family','rate','mg',timebins,'xlim',[-pretime duration]);
        [fit] = locfit(spikes,'alpha',[nSpikes/length(spikes) 0.15],'deg',2,'family','rate','mg',timebins,'xlim',[-pretime duration]);
        if doLog
            rate(j,:)=log(invlink(predict(fit,timebins,'band','n'),fit{4}{5})./(Ntrials-1))';
        elseif doSqrt
            rate(j,:)=sqrt(invlink(predict(fit,timebins,'band','n'),fit{4}{5})./(Ntrials-1))';
        else
            rate(j,:)=(invlink(predict(fit,timebins,'band','n'),fit{4}{5})./(Ntrials-1))';
        end            
    end; 
    close(wh);

    
    spikes = sort(cat(1,aligned{:}));
    spikes = [-pretime - 1; spikes; duration-pretime+1];
    [fit] = locfit(spikes,'alpha',[nSpikes/length(spikes) 0.15],'deg',2,'family','rate','mg',timebins,'xlim',[-pretime duration]);
    if doLog
        rateAll = log(invlink(predict(fit,timebins,'band','n'),fit{4}{5})./Ntrials)';
    elseif doSqrt
        rateAll = sqrt(invlink(predict(fit,timebins,'band','n'),fit{4}{5})./Ntrials)';
    else
        rateAll = (invlink(predict(fit,timebins,'band','n'),fit{4}{5})./Ntrials)';
    end

    ratem = mean(rate,1);
    ratese = sqrt((Ntrials-1).^2/Ntrials.*trimmean((rate - repmat(ratem,[Ntrials 1])).^2,5));
    
    rateAll = 2*rateAll - ratem; % modest bias correction
    rateci = [rateAll+tinv(0.975,Ntrials-1)*ratese; max(rateAll-tinv(0.975,Ntrials-1)*ratese,0)];

    if doLog
        ratem = exp(rateAll);
        rateci = exp(rateci);
        rate = exp(rate);
    elseif doSqrt
        ratem = (rateAll).^2;
        rateci = (rateci).^2;
        rate = (rate).^2;
    end
    
else % without the ability to jackknife, use a parametric estimate of the standard error
    spikes = cat(1,aligned{:});
    spikes = [-pretime - 1; spikes; duration-pretime+1];
    [fit] = locfit(spikes,'alpha',[nSpikes/length(spikes) 0.05*0.75/2],'deg',2,'family','rate','mg',timebins,'xlim',[-pretime duration]);
    ypp = predict(fit,timebins,'band','y');
    bands = ypp{3};
    lo = invlink(bands(:,1),fit{4}{5});
    up = invlink(bands(:,2),fit{4}{5});
    ratem = (invlink(ypp{1},fit{4}{5})./(Ntrials-1))';
    Ntrials = length(aligned);
    rateci = [up(:)'; lo(:)']/Ntrials;
    ratem = ratem/Ntrials;
    rate = ratem;
end;    
return;
