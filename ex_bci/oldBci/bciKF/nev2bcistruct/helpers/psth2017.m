function psth = psth2017(dat,alignendflag,endlimit)
    if nargin<2
        alignendflag = 0;
    end
    if nargin<3
        endlimit = zeros(length(dat),1);
    end
    
    
    lengthcount = zeros(length(dat),1);
    for n = 1:length(dat)
        lengthcount(n) = size(dat(n).counts,2);
    end
    channelcount = [{}];
    for n = 1:size(dat(n).counts,1)
        channelcount{n} = nan(length(dat),max(lengthcount));
    end
    
    for n = 1:length(dat)
        if endlimit(n) == 0
            spikecounts = dat(n).counts;
        else
            spikecounts = dat(n).counts(:,1:endlimit(n));
        end
        for m = 1:size(spikecounts,1)
            if alignendflag == 0
                channelcount{m}(n,1:size(spikecounts,2)) = spikecounts(m,:);
            else
                channelcount{m}(n,(end-size(spikecounts,2)+1):end) = spikecounts(m,:);
            end
        end
    end
    psth = zeros(size(spikecounts,1),max(lengthcount));
    for n = 1:length(channelcount)
        psth(n,:) = nanmean(channelcount{n},1);
    end
end