function spikenev = unpackSpikes(spikeinfo,spiketimesdiff,firstspike,Fs)
    newspiketimes = cumsum([double(firstspike);double(spiketimesdiff)])/Fs;
    spikenev = [double(spikeinfo) newspiketimes];
end