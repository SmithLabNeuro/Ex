function testJuice(drops,pams,dms)
% drops = number of trials 
% dms = duration of drops
% pams = pause between trials
counter = 0;
for n = 1:drops
    unixSendPulse(18,dms)
    pause(pams/1000)
counter = counter +1
end

