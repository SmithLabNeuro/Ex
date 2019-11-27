% This example will modulate the amplitude of a delivered stimulation pulse
% on electrode 1 of the first Micro+Stim connected based on input from a 
% sensor connected to analog input channel 1. In this example a force
% sensitive resistor is used where increasing force increases the voltage
% input on analog input channel 1. When the force level is below threshold,
% no stimulation will occur.  When the force level is above the threshold,
% the stimulation pulses will be delivered at a constant frequency with the
% amplitude of the biphasic pulse modulated based on the level of the 
% input.
%

%% Clean the world
close all; fclose('all'); clc; clear all; 


%% %%%%%%%%%%%%% Soft Parameter Initialization - May change as needed %%%%%%%%%%%%%%
anlgChan_idx  = 1;      % analog I/O FE channel to monitor for force
stimChan_idx  = 1:32;   % stim FE channel
stimFreq      = 10;     % stimulation Frequency (Hz)
phaseDur_us   = 100;    % Duration of cathodic and anodic phases of stim (must be mulitples of 33.3333333 for this eg.)
fs_us         = 200;    % duration of stim fast settle of recording (us)
window_ms     = 20;     % Length of window to analyze incoming analog data in control loop

%% %%%%%%%%%%%%%%%% Fixed Parameter Initialization - DO NOT CHANGE %%%%%%%%%%%%%%%%%
nipClock_us    = 1e6/3e4;           % 33.333 us
nipClock_ms    = nipClock_us * 1e3; % 0.0333 ms
msec2nip_clk   = 30;
nip_clk2sec    = 1/3e4;
nip_clk2msec   = nip_clk2sec * 1e3;
nip_clk2usec   = nip_clk2sec * 1e6;
AMP_NEURAL     = 0;              % used to set the input channel amp to measure neural voltages
AMP_STIM       = 1;              % used to set the input channel amp to measure stim voltage
stimAmp2V      = 0.50863e-3;
stimAmp2uV     = stimAmp2V * 1e6;

% Specific parameters used for this example
forceThresh_mV = 150;
forceMax_mV    = 4500;
forceRange_mV  = forceMax_mV - forceThresh_mV;

%% %%%%%%%%%%%%%%%%%% NIP Hardware Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize xippmex
status = xippmex;
if status ~= 1; error('Xippmex Did Not Initialize');  end

% Find all Micro+Stim and Analog channels
anlgChans = xippmex('elec','analog');
stimChans = xippmex('elec','stim');

% NOTE: This demo expects there to be one analog i/o and one micro+stim 
% front end attached to the NIP. So let's check if this is indeed the case.
if isempty(anlgChans); error('No analog hardware detected'); end

%make sure there is at least one micro+stim front end present
if isempty(stimChans); error('No stimulation hardware detected');  end

% Turn on 1 kS/s data stream, turn off 30 kS/s data stream to reduce load
xippmex('signal', anlgChans, '1ksps', ones(1,length(anlgChans)))
xippmex('signal', anlgChans, '30ksps', zeros(1,length(anlgChans)))

% Give the NIP some time to process any commands we have sent
pause(0.5)

%% %%%%%%%%%%%%% Stim Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Enable stimulation on the NIP.
% NOTE: If stimulation is not enabled, xippmex will enable and disable 
% stimulation for each stim sequence it receives. This will have the 
% undesirable side effect of killing queued stimulation sequences. Always 
% enable stimulation on the NIP before queueing multiple stimulation
% sequences.
xippmex('stim', 'enable', 0); pause(0.5)
xippmex('stim', 'enable', 1); pause(0.5)

% NOTE: This size is fixed and cannot be exceeded. Overflowing the buffer
% will cause the NIP to shut down all stimulation.
stimQueueSize = 8;

pw_cycls = floor(phaseDur_us / nipClock_us);
fs_cycls = floor(fs_us / nipClock_us);

for i = 1:32
    cmd(i).elec    = stimChans(stimChan_idx(i));
    cmd(i).period  = floor(1/nip_clk2sec / stimFreq);
    cmd(i).repeats = 200;
    cmd(i).action  = 'curcyc';  % Valid 'action' options are: ['immed', 'curcyc', 'allcyc', 'at-time']
    
    % Initialize cathodic phase of stim
    cmd(i).seq(1) = struct('length', pw_cycls, 'ampl', 0, 'pol', 0, ...
        'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', AMP_STIM);

    % Initialize interphase interval, i.e., time between cathodic an anodic 
    % phases of bipolar waveform. In this example, it is one clock cycle.
    cmd(i).seq(2) = struct('length', 1, 'ampl', 0, 'pol', 0, ...
        'fs', 0, 'enable', 0, 'delay', 0, 'ampSelect', AMP_STIM);

    % Initialize anodic phase of stim
    cmd(i).seq(3) = struct('length', pw_cycls, 'ampl', 0, 'pol', 1, ...
        'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', AMP_STIM);

    % Initialize fast settle after stimulation pulse
    if fs_cycls > 0
        cmd(i).seq(4) = struct('length', fs_cycls, 'ampl', 0, 'pol', 1, ...
            'fs', 1, 'enable', 1, 'delay', 0, 'ampSelect', AMP_STIM);
    end
    
    cmdClear(i).elec = stimChans(stimChan_idx(i));
    cmdClear(i).period  = 20;
    cmdClear(i).repeats = 1;
    cmdClear(i).action  = 'immed';

    cmdClear(i).seq(1) = struct('length', 3, 'ampl', 0, 'pol', 0, ...
        'fs', 0, 'enable', 0, 'delay', 0, 'ampSelect', AMP_NEURAL);

end
                     

%% %%%%%%%%%%%%%%%% Set control parameters for stimulation %%%%%%%%%%%
window_clk   = window_ms * msec2nip_clk;
lastNipTime  = 0;
nipOffTime   = 0;
stimOff      = 0;

% Frequency of control loop. In this example, this cannot be faster than
% stimulation frequency. This is to prevent buffer commands to 
% stack up for each electrode (>8), which would cause buffer overrun
% and shut down stimulation for that electrode.
loopPrd_s = 1/stimFreq;

%% %%%%%%%%%%%%%% Run Control Loop %%%%%%%%%%%%%%%%%%%%%%%%
loopTimer = tic; loopTimeCur = toc(loopTimer);
while true
    
    % wait until loop period has occured
    nextTime = loopTimeCur + loopPrd_s;
    while toc(loopTimer) < nextTime; end
    loopTimeCur = toc(loopTimer);
    
    % Get the current NIP time
    curNipTime = xippmex('time');
    
    % chek if the NIP is still on-line
    if curNipTime == lastNipTime
        if nipOffTime == 0
            tic;
        end
        nipOffTime = toc;
        % if the nip has been offline for a second abort the program
        if nipOffTime > 1
            xippmex('close');
            error('NIP appears to be off-line. Exiting program... Bye!')
        end
    % if the NIP is on-line clear the off timer
    else
        nipOffTime = 0;
    end
    % update NIP time for the next pass through the loop    
    lastNipTime = curNipTime;
    
    % collect recent analog data
    winStart_clk = curNipTime - window_clk;
    [anlgData, actTime] = xippmex('cont', anlgChans(anlgChan_idx), window_ms, '1ksps', winStart_clk);

    % if the force threshold has been crossed, modulate the amplitude of
    % stimulation.  if the force is below threshold, shut off stimulation
    meanAnlg = mean(anlgData);
    if meanAnlg > forceThresh_mV
        stimOff    = 0;
        forceDiff  = meanAnlg - forceThresh_mV;
        frqRngFrac = min(forceDiff/forceRange_mV, 1);
        for i = 1:32
            cmd(i).seq(1).ampl = floor (127 * frqRngFrac);        
            cmd(i).seq(3).ampl = floor (127 * frqRngFrac);
        end
        xippmex('stimseq', cmd);
                
    else
        if ~stimOff
            xippmex('stimseq', cmdClear);
            stimOff = 1;
        end
        continue
    end    
end

% let the user know the script is done
disp('Exiting Ripple Stim Demo .... Bye!')
xippmex('close');
         
         
         
         