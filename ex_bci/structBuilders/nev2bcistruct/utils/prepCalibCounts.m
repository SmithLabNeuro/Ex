function [countmat,trainingdat] = prepCalibCounts(trainingdat,starttimecode, endtimecode,bin)
countmat = [];
histchannels = [1:96 257:352];
for n = 1:length(trainingdat)
    % prep spikes and trial codes
    spikes = unpackSpikes(trainingdat(n).spikeinfo, trainingdat(n).spiketimesdiff,trainingdat(n).firstspike,30000);
    if length(unique(trainingdat(n).spikeinfo(:,2)))>1
        spikes = spikes(~ismember(spikes(:,2),[0 255]),:);
    end
    spikes(:,1) = spikes(:,1);
    trialcodes = trainingdat(n).trialcodes;
    starttime = trialcodes(trialcodes(:,2)==starttimecode,3);
    spiketimes = (spikes(:,3));
    spikechannels = spikes(:,1);
    endtime = trialcodes(trialcodes(:,2)==endtimecode,3);
    tempcounts = [];
    
    % compute spike counts and deal with boundary cases
    if ~isempty(starttime) && ~isempty(endtime)
        times = starttime(1):bin:endtime(1);
        if length(times)>1    
        tempcounts = histcn([spikechannels spiketimes],histchannels,times);
        if ~isempty(find(spiketimes==times(end),1))
            if size(tempcounts,2)>1
                tempcounts = tempcounts(:,1:(end-1));
            else
                tempcounts = [];
            end
        end
        if size(tempcounts,1)+1 == length(histchannels)
            tempcounts(end+1,:) = zeros(1,size(tempcounts,2)); 
        end
        trainingdat(n).counttimes = times;
        trainingdat(n).nBins = length(times)-1;
        end
    end

    % store counts
    trainingdat(n).counts = tempcounts;  
    countmat = [countmat tempcounts];
end
end