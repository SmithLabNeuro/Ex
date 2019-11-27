function spikeCounts = getSpikesCounts(spikeTimes,spikeChannels,binTimes,histChannels)
lbin = length(binTimes);
indices=zeros(lbin,1);
for n = 1:lbin
indices(n) = binaryTimeSearch(spikeTimes,binTimes(n));
end

spikeCounts = zeros(length(histChannels),lbin-1);
for n = 1:(lbin-1)
    spikeCounts(:,n)=histc(spikeChannels(indices(n):(indices(n+1)-1)),histChannels);
end
end