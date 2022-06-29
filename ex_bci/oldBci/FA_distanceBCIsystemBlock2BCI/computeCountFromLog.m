function inputdat = computeCountFromLog(inputdat)
nev_info = inputdat(1).nevinfo;
secondchannels = 128*2+1;
goodchannels = [1:96 secondchannels:secondchannels+95];
for n = 1:length(inputdat)
    if ~isempty(inputdat(n).binstarttimes)
        spikes = unpackSpikes(inputdat(n).spikeinfo, inputdat(n).spiketimesdiff,inputdat(n).firstspike,1);
        if length(unique(spikes(:,2)))>1
            spikes = spikes(~ismember(spikes(:,2),[0 255]),:);
            %spikes = spikes(~ismember(spikes(:,2),[255]),:);
        end
        niptime = nev_info.nevclockstart;
        spiketimes = spikes(:,3)+niptime;
        spikechans =  spikes(:,1);
        starttimes = [inputdat(n).binstarttimes; inputdat(n).binendtimes(end)];
        endtimes = inputdat(n).binendtimes;
        tempcounts = zeros(length(goodchannels),length(starttimes));
        %     for m = 1:length(goodchannels)
        %         for right = 1:size(tempcounts,2)
        %             temp=spiketimes(spikechans == goodchannels(m));
        %         tempcounts(m,right) = sum(temp>=starttimes(right)& temp <=endtimes(right) );
        %         end
        %     end
        tempcounts = histcn([spikechans spiketimes],goodchannels,starttimes);
        %tempcounts = getSpikesCounts(spiketimes,spikechannels,times,histchannels);
        if ~isempty(find(spiketimes==starttimes(end),1))
            if size(tempcounts,2)>1
                tempcountspre = tempcounts(:,1:(end-1));
                tempcountspre(:,end) = tempcountspre(:,end)+tempcounts(:,end);
                tempcounts = tempcountspre;
            else
                tempcounts = [];
            end
        end
        if size(tempcounts,1)+1 == length(goodchannels)
            tempcounts(end+1,:) = zeros(1,size(tempcounts,2));
        end
        inputdat(n).logcount = tempcounts;
    else
        inputdat(n).logcount = [];
    end
end
end
