%% test code for Rig 1

% addpath('C:\Documents and Settings\smithlab\My Documents\Dropbox\smithlab\Ex\xippmex');
% addpath('C:\Users\ripple\Dropbox\smithlab\Ex\xippmex');

status = xippmex;
if status~=1
    error('oh no');
end
% maybe check to see 1-32 are stim channels
stimChans  = xippmex('elec','stim')
% if (sum(diff(stimChans),1:32) > 0)
%     error('oh no');
% end
opers = xippmex('opers')

xippmex('signal',stimChans,'stim',stimChans)
%xippmex('signal',1:32,'stim',[1:32])

% build cmd - just prepares it
%stim_cmd = formatUstimCmd(1)
% [stim_cmd] = xippmexStimCmd(stim_chan,pulse_width,stim_freq,stim_dur,stim_amp)
stim_chan = stimChans(1);
pulse_width = 250;
stim_freq = 350;
stim_dur = 200;
stim_amp = 100;

[stim_cmd] = xippmexStimCmd(stim_chan,pulse_width,stim_freq,stim_dur,stim_amp)


%% send stim cmd - actually stimulates!
xippmex('stimseq',stim_cmd)

%% figure out the current step setting
xippmex('stim','res',stim_chan)

%% read in some spikes
status = xippmex;
[count, timestamps, waveforms, units] = xippmex('spike', [1], 0);

clear count timestamps waveforms units y;
for I=1:10
    tic;
    [count, timestamps, waveforms, units] = xippmex('spike', [1:32], 0);
    y(I)=toc;
end
disp(mean(y));

%%% MATT & ROMA STOPPED HERE

%% call xippmex once to deal with initialization lag
if isempty(TrialCount)
    TrialCount = 1;
    xippmex
else
    TrialCount = TrialCount + 1;
end




% set up microstim parameters
[stim_cmd] = GenerateUstimCommands;



%send ustim command at right part of ex_ code
                   xippmex('stimseq', stim_cmd);
                   
                
% try to crash xippmex
matlabUDP('close'); %first close any open UDP sessions. -ACS 13AUG2013
matlabUDP('open','192.168.1.11','192.168.1.10', 4243)

% begin the trial

status = xippmex;
stimChans=xippmex('elec','stim');
StimAmp = [50 50 50];
StimChan = [1 2 3];
if ~isscalar(StimAmp), StimAmp(StimChan<1) = []; end;
StimChan(StimChan<1) = []; 
if any(setdiff(StimChan,stimChans))
    error('Unable to stimulate on requested channel %d',setdiff(StimChan,stimChans));
end
if ~isempty(StimChan) && any(StimAmp),
    xippmex('signal',StimChan,'stim',StimChan); 
    stim_cmd = xippmexStimCmd(StimChan,250,350,200,StimAmp);
end;

xippmex('stimseq',stim_cmd);
pause(1);

% end of trial

%xippmex('close');
matlabUDP('close');

pause(0.2); disp('between runs');

matlabUDP('close'); %first close any open UDP sessions. -ACS 13AUG2013
matlabUDP('open','192.168.1.11','192.168.1.10', 4243)


%status = xippmex;
stimChans=xippmex('elec','stim');
StimAmp = [50 50 50];
StimChan = [1 2 3];
if ~isscalar(StimAmp), StimAmp(StimChan<1) = []; end;
StimChan(StimChan<1) = []; 
if any(setdiff(StimChan,stimChans))
    error('Unable to stimulate on requested channel %d',setdiff(StimChan,stimChans));
end
if ~isempty(StimChan) && any(StimAmp),
    xippmex('signal',StimChan,'stim',StimChan); 
    stim_cmd = xippmexStimCmd(StimChan,250,350,200,StimAmp);
end;

xippmex('stimseq',stim_cmd);
pause(1);
%xippmex('close');
matlabUDP('close');

                   
%%%%%%%%%%%%%%%%
% smallest test possible
matlabUDP('close'); %first close any open UDP sessions. -ACS 13AUG2013
matlabUDP('open','192.168.1.11','192.168.1.10', 4243)

StimChan=xippmex('elec','stim');
StimAmp = [50 * ones(length(StimChan),1)]';
stim_cmd = xippmexStimCmd(StimChan,250,350,50,StimAmp);

status = xippmex

for I=1:10
    %status = xippmex
    StimChan=xippmex('elec','stim');
    xippmex('signal',StimChan,'stim',StimChan);
    xippmex('stimseq',stim_cmd);
    pause(0.1);
    %xippmex('close')
end

matlabUDP('close'); %first close any open UDP sessions. -ACS 13AUG2013
