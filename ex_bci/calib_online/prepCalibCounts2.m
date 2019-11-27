function [countmat, anglemat, distancemat,trainingdat,Xmat,Ymat] = prepCalibCounts(trainingdat,starttimecode, endtimecode,minchannelnum,bin,maxbins)
if ~exist('maxbins','var')
    maxbins = 100000; % huge number that should never be exceeded
end


angletemp = driftchoiceextractparam(trainingdat,'angle=');
unangle = unique(angletemp);
% if unangle(1) ~= 0
%     angletemp = angletemp - unangle(1); % subtraction needed for some data sets;
% end
distancetemp = driftchoiceextractparam(trainingdat,'distance=');
angle = angletemp;%(distancetemp>80);
distance = distancetemp;%130*ones(length(distancetemp),1);


countmat = [];
anglemat = [];
Xmat = [];
Ymat = [];
distancemat = [];
for n = 1:length(trainingdat)
    spikes = trainingdat(n).spikes(trainingdat(n).spikes(:,1)>minchannelnum,:);
    spikes(:,1) = spikes(:,1)-minchannelnum;
    trialcodes = trainingdat(n).trialcodes;
    starttime = trialcodes(trialcodes(:,2)==starttimecode,3);
    
    endtime = trialcodes(trialcodes(:,2)==endtimecode,3);
    %times = [starttime (starttime+bin)];
    %times = [(endtime-bin) endtime];
    times = starttime(1):bin:endtime(1);
    tempcounts = [];
    for m = 1:(length(times)-1)
        if m > maxbins
            break
        end
        temp = histc(spikes(spikes(:,3)>=times(m) & spikes(:,3)<times(m+1),1),[1:96 257:352]);
        if size(temp,1)<size(temp,2)
            temp = temp';
        end
        tempcounts = [tempcounts temp];%20*temp/norm(temp)];%
    end
    trainingdat(n).counts = tempcounts;
    
    theta = deg2rad(angle(n));
    trainingdat(n).newX = round(distance(n)*cos(theta));
    trainingdat(n).newY = round(distance(n)*sin(theta));
    trainingdat(n).angle = angle(n);
    countmat = [countmat tempcounts];
    anglemat = [anglemat angle(n)*ones(1,size(tempcounts,2))];
    Xmat = [Xmat trainingdat(n).newX*ones(1,size(tempcounts,2))];
    Ymat = [Ymat trainingdat(n).newY*ones(1,size(tempcounts,2))];
    distancemat = [distancemat distance(n)*ones(1,size(tempcounts,2))];
    trainingdat(n).counttimes = times;
    trainingdat(n).nBins = length(times)-1;
end
end