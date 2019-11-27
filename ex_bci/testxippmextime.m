status = xippmex;
elec = xippmex('elec','nano');
numloops = 1000;
looptime = zeros(numloops,1);
spikecounts = zeros(numloops,1);
[spkCount, spkTimestamps, spkWaveforms] = xippmex('spike', elec,0); % get spikes
for n = 1:numloops
    tic;
    timenew=xippmex('time'); % set time for loop n+1
    [spkCount, spkTimestamps, spkWaveforms] = xippmex('spike', elec,0); % get spikes
    spikecounts(n) = sum(spkCount);
    % ==================
    % analysis code here    
    % ==================

    t=xippmex('time');
    while t-timenew < 0.0001*30000
        t=xippmex('time');
    end
    looptime(n) = toc;
end
looptime
corr(looptime, spikecounts)
xippmex('close')
figure
scatter(looptime,spikecounts)