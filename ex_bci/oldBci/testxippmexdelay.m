% This script estimates the minimum possible spike time on each channel
% (+ a small constant lag) and plots the offset from the minimum spike time on
% channel 1. The plot shows a delay of minimum spike time that
% increases as channel number increases.
clear all
%addpath(genpath('../xippmex-1.3.1.64_linux'));
%addpath('/home/smithlab/Dropbox/smithlabrig/Ex/xippmex-1.3.1.64_linux')
addpath('/home/smithlab/Dropbox/smithlabrig/Ex/xippmex-1.11')
status = xippmex; % start xippex
elec = xippmex('elec','nano'); %get valid electrodes
numloops = 1000; % this needs to be large if spike rate is low in order to estimate min spike time on each channel
vals = nan(length(elec),numloops); % use nan to ignore empty channels/bad spikes
time=xippmex('time');%set first time
xippmex('spike', elec);% clear xippmex buffer
for n = 1:numloops;
    tic;
    if mod(n,100)==0 %print out progress for my sanity
    fprintf('Loop %i\n',n);
    end
    timenew=xippmex('time'); % set time for loop n+1
    [spkCount, spkTimestamps, spkWaveforms] = xippmex('spike', elec); % get spikes
    empty=cellfun(@(x)(isempty(x)),spkTimestamps); % find channels with no spikes
    nonemptytimestamps = spkTimestamps(~empty); % remove empty channels
    temp1=(cellfun(@(x)(min(x)),nonemptytimestamps)); % find min spike time on each channel
    if sum(temp1==0)>0
        display('zero val')
    end
    temp1(temp1==0)=nan; % remove cases where spike time is zero (why does this happen?)
    vals(~empty,n) = (temp1-time)/30000; % subtract off fixed time and divide by sampling rate
    time = timenew; % 
    while toc<0.05
    end
end
temp=nanmin(vals,[],2);
figure; % plot minimum times from offset
plot(temp-temp(1));%ylim([-0.001 0.002]);
xlabel('Channel Number');ylabel('Difference from Channel 1 Min Time')
xippmex('close'); % close xippmex