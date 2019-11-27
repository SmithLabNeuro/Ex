function dat = getNS5Data(dat,filepath,ns5filename,filetype)
%% figure out ns5 start times
if iscell(ns5filename)
    starttimes=zeros(length(dat),1);for n =1:length(dat);starttimes(n)=dat(n).nevinfo.nevclockstart;end
    unstarttimes = unique(starttimes);
    if length(unstarttimes)~= length(ns5filename)
        error('Must have same number of NS5 files as NEV files')
    end
    for n = 1:length(ns5filename)
        data(n) = read_nsx([filepath,ns5filename{n},filetype],'chanindx',1:4);
    end
else
     data = read_nsx([filepath,ns5filename,filetype],'chanindx',1:4);
     unstarttimes = dat(1).nevinfo.nevclockstart;
end
 dataFs = double(data(1).hdr.Fs);
 %clockFs = double(data.hdr.clockFs);
 downsampleeye = 30;
 downsamplediode = 3;
 for tind = 1:length(dat)
     thisstart = dat(tind).nevinfo.nevclockstart;
     dataind = find(unstarttimes==thisstart);
     msec = dat(tind).trialcodes(:,3);
     codes = dat(tind).trialcodes(:,2);
     %codesamples = nevMsec2nsxSample(msec,dataFs,data.hdr.timeStamps,clockFs);
     codesamples = round(msec*dataFs);
     eyedata.codesamples = [codes codesamples];
     eyedata.trial = downsample(data(dataind).data(1:2,codesamples(codes==1):codesamples(codes==255))',downsampleeye)';
     eyedata.startsample = codesamples(1)/downsampleeye;
     eyedata.dataFs = dataFs/downsampleeye;
     diode.codesamples = [codes codesamples];
     diode.trial = int16(downsample(data(dataind).data(3,codesamples(1):codesamples(end)),downsamplediode));
     diode.startsample = codesamples(1)/downsamplediode;
     diode.dataFs = dataFs/downsamplediode;
     
     pupil.codesamples = [codes codesamples];
     pupil.trial = int16(downsample(data(dataind).data(4,codesamples(1):codesamples(end)),downsamplediode));
     pupil.startsample = codesamples(1)/downsamplediode;
     pupil.dataFs = dataFs/downsamplediode;
     
     dat(tind).eyedata = eyedata;
     dat(tind).diode = diode;
     dat(tind).pupil = pupil;
end
end