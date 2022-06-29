function [countmat,trainingdat] = prepCalibCounts(trainingdat,starttimecode, endtimecode,minchannelnum,bin,maxbins)
    if ~exist('maxbins','var')
        maxbins = 100000; % huge number that should never be exceeded
    end

%% old way of getting angle and distance
% angletemp = driftchoiceextractparam(trainingdat,'angle=',2); % set to 2 b/c angle is changed mid-trial in some tasks
% unangle = unique(angletemp);
% % if unangle(1) ~= 0
% %     angletemp = angletemp - unangle(1); % subtraction needed for some data sets;
% % end
% distancetemp = driftchoiceextractparam(trainingdat,'distance=');
% angle = angletemp;%(distancetemp>80);
% distance = distancetemp;%130*ones(length(distancetemp),1);

    %%

    countmat = [];
    histchannels = [1:96 257:352];
    for n = 1:length(trainingdat)
        %goodspikes = trainingdat(n).spikeinfo(:,1)>minchannelnum;
        spikes = unpackSpikes(trainingdat(n).spikeinfo, trainingdat(n).spiketimesdiff,trainingdat(n).firstspike,30000);
        if length(unique(trainingdat(n).spikeinfo(:,2)))>1
            spikes = spikes(~ismember(spikes(:,2),[0 255]),:);
            %spikes = spikes(~ismember(spikes(:,2),[255]),:);
        end
        %[double(trainingdat(n).spikeinfo) double(trainingdat(n).spiketimes(goodspikes))];
        spikes(:,1) = spikes(:,1)-minchannelnum;
        trialcodes = trainingdat(n).trialcodes;
        starttime = trialcodes(trialcodes(:,2)==starttimecode,3);
        spiketimes = (spikes(:,3));
        spikechannels = spikes(:,1);
        endtime = trialcodes(trialcodes(:,2)==endtimecode,3);
        %times = [starttime (starttime+bin)];
        %times = [(endtime-bin) endtime];
        tempcounts = [];
        if ~isempty(starttime) && ~isempty(endtime)
            times = starttime(1):bin:endtime(1);
            if length(times)>1
    %         tempcounts = zeros(length(histchannels),length(times));
    %         for k = 1:length(histchannels)
    %             temptimes = spiketimes(spikechannels == histchannels(k));
    %             tempcounts(k,:) = histc(temptimes,times);
    %         end
    %         
            tempcounts = histcn([spikechannels spiketimes],histchannels,times);
%             tempcounts = getSpikesCounts(spiketimes,spikechannels,times,histchannels);
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
    %         if size(tempcounts,2)~=length(times)
    %             display('counts doing weird stuff')
    %         end
    %         tempcounts = tempcounts(:,end-1);
            end
        end

        trainingdat(n).counts = tempcounts;
        countmat = [countmat tempcounts];

    end
end