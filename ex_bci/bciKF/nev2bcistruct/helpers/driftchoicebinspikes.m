function dat = driftchoicebinspikes(dat,event1,event2,binsize,channelinds,offsetval)
if nargin >5
    offsetval = offsetval;
else
    offsetval = 0;
end

for n =1:length(dat)
    if mod(n,100)==0; fprintf('Binning counts for trial %d of %d...\n',n,length(dat));end
    bintimesstart = dat(n).events(dat(n).events(:,2)==event1,3);
    bintimesend = dat(n).events(dat(n).events(:,2)==event2,3);
    if offsetval ~= 0
        bintimes = bintimesstart:offsetval:bintimesend;%-binsize;
        filtb = ones(1,binsize/offsetval);
    else
        bintimes = bintimesstart:binsize:bintimesend;
    end
    if ~isempty(bintimes)
        for m = 1:length(channelinds)
            if offsetval ~= 0
                tempbins = histc(dat(n).spikes(dat(n).spikes(:,1)==channelinds(m),3),bintimes);
                filttempbins = filter(filtb,1,tempbins(1:end-1));
                dat(n).spikecounts(m,:) = filttempbins(length(filtb):end)';
            else
                tempbins = histc(dat(n).spikes(dat(n).spikes(:,1)==channelinds(m),3),bintimes);
                dat(n).spikecounts(m,:) = tempbins(1:end-1);
            end
        end  
    else
        dat(n).spikecounts = [];
    end
end

end