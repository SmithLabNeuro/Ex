function [smoothedbins,startsteps,endsteps] = smoothSparse(timestamps,indata,nevflag,stepsizemin, window,mintime,maxtime)
timetoaverage = 1;%seconds Note that these values depending on the units of time in timestamps
stepsize = 1;% second
%stepsizemin = 60 * 1; % time in seconds
%window = 60 * 20; % time in seconds
    checktimesend = mintime:stepsize:maxtime;
    checktimesstart = checktimesend-timetoaverage;
if nevflag == 1
    channels = unique(indata(:,1));
    spiketimes = indata(:,3);
    spikechannels = indata(:,1);
    spikecountbins = histcn([spikechannels spiketimes],channels,[checktimesstart checktimesend(end)]);
    %spikecountbins = getSpikesCounts(spiketimes,spikechannels,[checktimesstart checktimesend],channels);
else
    spikecountbins = nan(1,length(checktimesend));
    for n = 1:length(checktimesend)
        spikecountbins(n) = mean(mean(indata(checktimesstart(n)<timestamps &checktimesend(n)>=timestamps)));
    end       

end

if nevflag == 1
    sumoverchannel = sum(spikecountbins,1);
    spikecountbinsnan = spikecountbins;
    spikecountbinsnan(:,sumoverchannel==0)=nan(size(spikecountbins,1),sum(sumoverchannel==0));
    endsteps = (window+1):stepsizemin:size(spikecountbins,2);
    startsteps = endsteps-window;
    smoothedbins = zeros(size(spikecountbins,1),length(endsteps));
else
    endsteps = (window+1):stepsizemin:length(spikecountbins);
    startsteps = endsteps-window;
    smoothedbins = zeros(1,length(endsteps));
end
for n = 1:length(endsteps)
    if nevflag == 1
        tempbin = nanmean(spikecountbinsnan(:,startsteps(n):endsteps(n)),2);
        if isempty(find(isnan(tempbin),1))
            smoothedbins(:,n) = tempbin;
        end
    else
        tempbin = nanmean(spikecountbins(startsteps(n):endsteps(n)));
        if isempty(find(isnan(tempbin),1))
            smoothedbins(n) = tempbin;
        end
    end
end
end
