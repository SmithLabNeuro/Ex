function [ bciDat ] = getOnlineCounts( bciDat )
%GETBCICOUNTS Summary of this function goes here
%   Detailed explanation goes here

    rmTrials = false(length(bciDat),1);
    for n = 1:length(bciDat)
        % get spiking info
        firstspike = bciDat(n).firstspike + bciDat(n).nevinfo.nevclockstart;
        spiketimes = cumsum([firstspike; double(bciDat(n).spiketimesdiff)]);
        spikeinfo = double(bciDat(n).spikeinfo(:,1));
        
        % get sort info
        sortCodes = double(bciDat(n).spikeinfo(:,2));
        sortIdx = ~ismember(sortCodes,[0 255]);
        spiketimes_sorted = spiketimes(sortIdx);
        spikeinfo_sorted = spikeinfo(sortIdx);
        
        % throw out trials where the bci didn't start
        if isempty(bciDat(n).binstarttimes)
            rmTrials(n) = true;
            continue;
        end
        
        % compute counts
        binEdges = [bciDat(n).binstarttimes; bciDat(n).binendtimes(end)];
        bciDat(n).T = length(binEdges)-1;
        channelNums = [1:96 257:352];
        counts = histcn([spiketimes spikeinfo],binEdges,channelNums)';
        if size(counts,2)>(length(binEdges)-1)
            counts(:,end-1) = counts(:,end-1)+counts(:,end);
            counts(:,end) = [];
        end
        bciDat(n).onlineCounts = counts;
        
        % compute sorted counts
        sortedCounts = histcn([spiketimes_sorted spikeinfo_sorted],binEdges,channelNums)';
        if size(sortedCounts,2)>(length(binEdges)-1)
            sortedCounts(:,end-1) = sortedCounts(:,end-1)+sortedCounts(:,end);
            sortedCounts(:,end) = [];
        end
        bciDat(n).onlineSortedCounts = sortedCounts;
    end
    bciDat = bciDat(~rmTrials);

end

