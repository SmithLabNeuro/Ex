function [goodcounts, percentexclude] = excludeCrossChannelNoise(spikecell,timebin,numevents)

if nargin==1
    numevents = 5;
    timebin = 0.001;
end

%% convert spike cell to spikecodes
channelnum=[];
spiketimes = [];
for n = 1:length(spikecell)
    if isrow(spikecell{n})
        spiketimestemp = spikecell{n}';
    else
        spiketimestemp = spikecell{n};
    end
    channelnumtemp = n*ones(length(spiketimestemp),1);
    spiketimes = [spiketimes;spiketimestemp];
    channelnum = [channelnum;channelnumtemp];
end

[spiketimes, sortind] = sort(spiketimes,'ascend');
channelnum = channelnum(sortind);
%% find spikes to exclude

diftime = [1;spiketimes(2:end)-spiketimes(1:end-1)]; % augment with 1 so first spike isn't dropped
filteredcoinpre = filter(ones(1,numevents),1,diftime);
filtcoincidences = filteredcoinpre<timebin;
excludecell = ones(length(spiketimes),1);
inds = find(filtcoincidences);
for n = 1:length(inds)
    ind1 = inds(n)-numevents+1;
    if ind1>0
    excludecell(ind1:inds(n)) = 0;
    end
end

goodcounts = zeros(length(spikecell),1);
for n = 1:length(goodcounts)
    goodcounts(n) = sum(excludecell(channelnum==n));
end
percentexclude = 1-sum(excludecell)/length(excludecell);
end