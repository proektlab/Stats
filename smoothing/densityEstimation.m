function [densitym,densityci,timebins,density] = densityEstimation(samples,limits,binsize,Nevents);
%function [densitym,densityci,timebins,density] = densityEstimation(samples,limits,binsize,Nevents);
%
% INPUTS:
%     samples, a cell array of spike times, one cell per trial
%     duration, a scalar, specifying the total amount of time in each entry of
%           the samples structure in seconds
%     pretime, a scalar, indicating the amount of "negative time" included
%           in the samples structure, in seconds
%     binsize, a scalar, specifying the time in seconds between the
%           interpolated density values
%     Nevents, a parameter that governs the bandwidth of the smoothing; the
%           nearest neighbor bandwidth will include this many events
%
% OUTPUTS:
%     densitym, the mean Poisson firing density
%     densityci, the jackknife confidence intervals for the Poisson firing
%           density
%     timepts, the time grid corresponding to each entry in densitym
%     density, the jackknife population of density estimates
if nargin<2, limits = []; end;
if nargin<3, binsize = []; end;
if nargin<4, Nevents = []; end;

if isempty(limits), limits =  [min(cat(1,samples{:})) max(cat(1,samples{:}))]; end;
if isempty(binsize), binsize = 0.001; end;
if isempty(Nevents), Nevents = 100; end;

if length(limits) == 1
    limits = [0 limits];
end;
limits = sort(limits);

doLog = 0;
doSqrt = 1; % note: sqrt is the variance-stabilizing link for a poisson distribution, and hence for densities, which are counts
if iscell(samples)
    if length(cat(1,samples{:}))>2
        Nsamples = length(samples);
    else
        Nsamples = 0;
    end
else
    Nsamples = 1;
    a{1} = samples;
    samples = a;
    clear a;
end

duration = limits(2)-limits(1);
timebins = (binsize.*(1:round(duration./binsize))+limits(1))';
if Nsamples>3
    wh=waitbar(0,'Local Density Estimation');
    for j=1:Nsamples 
        waitbar(j/Nsamples,wh);
        events = sort(cat(1,samples{[1:j-1 j+1:Nsamples]}));
        if ~isempty(events)
            [fit] = locfit(events,'alpha',[0.2 0.6],'deg',2,'family','qrdensity','mg',timebins,'xlim',limits);
            % [fit] = locfit(events,'alpha',[Nevents/length(events) 0.075],'deg',2,'family','qrdensity','mg',timebins,'xlim',limits);
            if doLog
                density(j,:)=log(invlink(predict(fit,timebins,'band','n'),fit{4}{5})./(Nsamples-1))';
            elseif doSqrt
                density(j,:)=sqrt(invlink(predict(fit,timebins,'band','n'),fit{4}{5})./(Nsamples-1))';
            else
                density(j,:)=(invlink(predict(fit,timebins,'band','n'),fit{4}{5})./(Nsamples-1))';
            end
        else
            density = NaN*ones(Nsamples,length(timebins));
        end
    end; 
    close(wh);
    densitym = mean(density,1);
    densityse = sqrt((Nsamples-1)/Nsamples.*sum((density - repmat(densitym,[Nsamples 1])).^2,1));
    densityci = [densitym+tinv(0.975,Nsamples-1)*densityse; max(densitym-tinv(0.975,Nsamples-1)*densityse,0)];
    if doLog
        densitym = exp(densitym);
        densityci = exp(densityci);
        density = exp(density);
    elseif doSqrt
        densitym = (densitym).^2;
        densityci = (densityci).^2;
        density = (density).^2;
    end
else % without the ability to jackknife, use a parametric estimate of the standard error
    events = cat(1,samples{:});
    if ~isempty(events)
        [fit] = locfit(events,'alpha',[Nevents/length(events) 0.05*0.75/2],'deg',2,'family','density','mg',timebins,'xlim',[limits(1) duration]);
        ypp = predict(fit,timebins,'band','y');
        bands = ypp{3};
        lo = invlink(bands(:,1),fit{4}{5});
        up = invlink(bands(:,2),fit{4}{5});
        densitym = (invlink(ypp{1},fit{4}{5})./(Nsamples))';
        Nsamples = length(samples);
        densityci = [up; lo]/Nsamples;
        densitym = densitym/Nsamples;
        density = densitym;
    else
        densitym = NaN.*timebins';
        densityci = NaN.*[timebins timebins]';
        density = densitym;
    end
end;    
return;
