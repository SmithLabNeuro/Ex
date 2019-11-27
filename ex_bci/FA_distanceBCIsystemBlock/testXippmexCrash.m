addpath(genpath('../../xippmex-1.11'))
microornano= 'nano';
status = xippmex;
elecs = xippmex('elec',microornano);
count = 0;
while 1
    %pause(0.001)
    count = count +1;
    aligntime = xippmex('time');
    [countbuffertemp,spiketimes,waveforms,units]=xippmex('spike',elecs,[zeros(1,length(elecs))]);
    if mod(count,1000)==0
        display(count)
    end
end