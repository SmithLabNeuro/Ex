
status = xippmex;
[countbuffertemp,spiketimes,~,units]=xippmex('spike',elecs,[zeros(1,length(elecs))]);
timestamp = min([spiketimes{:}]);
while 1
    pause(0.05)
    [countbuffertemp,spiketimes,~,units]=xippmex('spike',elecs,[zeros(1,length(elecs))]);
    minspiketime = min([spiketimes{:}]);
    maxspiketime = max([spiketimes{:}]);
    if minspiketime < timestamp || minspiketime==0
        fprintf('minspiketime=%d, maxspiketime=%d, timestamp=%d\n',minspiketime,maxspiketime, timestamp)
    end
    timestamp = minspiketime;
end
xippmex('close');
