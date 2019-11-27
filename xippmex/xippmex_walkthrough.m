% This is a Xippmex test script that will run through the following
% operations:
%
% 1) Initialize Xippmex
% 2) Find all Micro+Stim, Micro, Micro HV, and Nano channels (if they exist)
% 3) Configure the data streams on the Front Ends
% 4) Set filters for the Front End LFP and HiRes data streams
% 5) Configure Trellis information if Trellis in use on 2nd Computer 
%        %%%% NOTE: In Trellis, File Save 'Enable Remote Control' must be selected %%%%
% 6) Create and execute a stimulation pattern for a Micro+Stim
% 7) Collect data
% 8) Plot all collected data included delivered stimulation
% 9) Close Xippmex


%% Initializations 
% Clean the world
close all; fclose('all'); clc; clear all; 

% Initialize xippmex
status = xippmex;
if status ~= 1; error('Xippmex Did Not Initialize');  end

% Find all Micro+Stim and Micro/Nano channels
stimChans  = xippmex('elec','stim');
recChans   = [xippmex('elec','micro') xippmex('elec','nano')];

% Get NIP clock time right before turning streams on (30 kHz sampling)
timeZero = xippmex('time');

%% Data Stream Activation
% Configure electrode 5&6 on the Micro+Stim to have Stim data stream 
% active, all Hi-Res streams to active, and deactivate all other streams
if ~isempty(stimChans)
    xippmex('signal',stimChans(5:6),        'stim',   [1 1]);
    xippmex('signal',stimChans([1:4 7:32]), 'stim',   zeros(1,30));
    xippmex('signal',stimChans,             'hi-res', ones(1,length(stimChans)));   
    xippmex('signal',stimChans,             'raw',    zeros(1,length(stimChans)));
    xippmex('signal',stimChans,             'lfp',    zeros(1,length(stimChans)));
    xippmex('signal',stimChans,             'spk',    zeros(1,length(stimChans)));
end

% Configure the first 8 recording channels to have Spike streams active
% LFP active on all, and deactivate all other streams
if ~isempty(recChans)
    xippmex('signal',recChans(1:8),  'spk',    ones(1,8));
    xippmex('signal',recChans(9:32), 'spk',    zeros(1,24)); 
    xippmex('signal',recChans,       'lfp',    ones(1,length(recChans)));   
    xippmex('signal',recChans,       'raw',    zeros(1,length(recChans)));
    xippmex('signal',recChans,       'hi-res', zeros(1,length(recChans)));
end

%% Filter Settings
% Set the Hi-Res notch filter for the first Micro+Stim FE to use
% the 60/120/180 comb filter
if ~isempty(stimChans)
    [selection, filters] = xippmex('filter','list',stimChans(1),'hires notch');
    filterNum = find(strcmp({filters.label}, '60/120/180 Hz'));
    xippmex('filter','set', stimChans(1),'hires notch', filterNum);
end

% Set the LFP filter for the first recording FE to use the wide-band 0.3-250 Hz filter
if ~isempty(recChans)
    [selection, filters] = xippmex('filter','list',recChans(1),'lfp');
    filterNum = find(strcmp({filters.label}, '0.3-250 Hz'));
    xippmex('filter','set', recChans(1),'lfp', filterNum);
end

%% If a Trellis instance was found, start recording for a short segment
% NOTE: In Trellis, File Save 'Enable Remote Control' must be selected
% Get running Trellis instances. If one is running get the trial
% descriptor, set it to run for 10 s, and then start recording.
opers = xippmex('opers');
if ~isempty(opers)
    desc = xippmex('trial', opers(1));
    xippmex('trial', opers(1), 'stopped', desc.filebase, 10, 1);
    xippmex('trial', opers(1), 'recording');
end

%% Stimulation Settings
% Create a stimulation pattern for electrode 5&6 on the Micro+Stim with 
% bipolar stimulation between the electrodes, with a frequency of 30 Hz and
% a train length of two seconds. In this case, electrode 5 will have the 
% cathodic phase first and electrode 10 will have the anodic phase first. 
% So for the first phase, current will flow from electrode 5 into 
% electrode 6. 
if ~isempty(stimChans)  
    % stimulation parameters
    stimElectrodes         = stimChans(5:6);
    trainLength_ms         = [2000 2000];
    frequency_Hz           = [30 30];
    phaseDuration_ms       = [0.4 0.4];
    phaseAmplitude_steps   = [20 20];
    electrodeDelay_ms      = [0 0];
    polarity               = [1 0];
     
    % Generate stimulation string
    stimString = stim_param_to_string(stimElectrodes, trainLength_ms, frequency_Hz, ...
        phaseDuration_ms, phaseAmplitude_steps, electrodeDelay_ms, polarity);
    
    % Execute stimulation
    xippmex('stim',stimString);
end


%% Wait for recording to end on Trellis side if Trellis operators exist, otherwise pause for 5 s 
if ~isempty(opers)
    desc = xippmex('trial', opers(1));
    while ~strcmp(desc.status, 'stopped')
        pause(1); desc = xippmex('trial', opers(1));
    end
else
    pause(5)
end

%% Collect Data
if ~isempty(stimChans)
    [stimCount, stimTimestamps, stimWaveforms] = xippmex('spike', stimChans(5:6), [1 1]);        % stim waveforms
    [HiResData, HiResTimestamp]                = xippmex('cont',stimChans(9:12),5000,'hi-res');  % hires data
end

if ~isempty(recChans)
    [spkCount, spkTimestamps, spkWaveforms] = xippmex('spike', recChans(1:8), zeros(1,8)); % spike waveforms on Micro
    [LFPData, LFPTimestamp]                 = xippmex('cont',recChans(1:8),5000,'lfp');    % lfp data
end

%% Plot Stim Waveforms
if ~isempty(stimChans)
    stims = zeros(3e4*2.5,2);                          % create array of zeros to append stim waveforms to 
    startTime = stimTimestamps{1,1}(1) - timeZero;     % offset timestamps based on initial NIP clock sample
    
    % loop through each channel, then loop through each stimulation event
    % noted by a stimCount/stimTimestamp. Find when each stim event occured
    % in relation to the first stim event (stimWindow) and append the
    % extracted stimWaveform to the stims array.
    for i = 1:2        
        for j = 1:stimCount(i)
            stimWindow = (1:52) + (stimTimestamps{i,1}(j) - stimTimestamps{1,1}(1));         
            stims(stimWindow,i) = stims(stimWindow,i) + (cell2mat(stimWaveforms{i,1}(j)))';
        end      
    end
    
    time = (startTime:startTime+size(stims,1)-1)/3e4;  % create array of time for plot usage
    figure(); plot(time,stims(:,1)/1e3, time,stims(:,2)/1e3);
    xlabel('Time (s)'); ylabel('Stim Voltage (mV)'); title('Stimulation Waveforms'); legend('Stim Chan 5','Stim Chan 6'); 
end

%% Plot HiRes Data
if ~isempty(stimChans)    
    startTime = floor((HiResTimestamp - timeZero) / 3e4 * 2e3);   % determine when the first HiRes data sample was in relation to the first NIP clock time pulled earlier.
    time = (startTime:startTime+5*2e3-1)/2e3;                     % create array of time for plot usage
    figure(); plot(time,HiResData(1,:), time,HiResData(2,:), time,HiResData(3,:), time,HiResData(4,:));
    xlabel('Time (s)'); ylabel('Recorded Voltage (uV)'); title('Hi-Res Data'); legend('Stim Ch 9','Stim Ch 10','Stim Ch 11','Stim Ch 12');    
end

%% Plot Spike Data for Channel 3 on the Micro, HV, or Nano (if there are any spikes)
if ~isempty(recChans)
    time = (0:51)/30;  % 52 sample waveform at 30kHz sampling
    if spkCount(3)
        figure(); hold on; 
        for i = 1:spkCount(3); plot(time,cell2mat(spkWaveforms{3}(i))); end
        hold off; xlabel('Time (ms)'); ylabel('Recorded Voltage (uV)'); title('Spike Waveforms for Rec Ch 3');
    end
end

%% Plot LFP Data
if ~isempty(recChans)
    startTime = floor((LFPTimestamp - timeZero) / 3e4 * 1e3);  % determine when the first HiRes data sample was in relation to the first NIP clock time pulled earlier.
    time = (startTime:startTime+5*1e3-1)/1e3;                  % create array of time for plot usage
    figure(); plot(time,LFPData(1,:), time,LFPData(2,:), time,LFPData(3,:), time,LFPData(4,:));
    xlabel('Time (s)'); ylabel('Recorded Voltage (uV)'); title('LFP Data'); legend('Rec Ch 1','Rec Ch 2','Rec Ch 3','Rec Ch 4');    
end


%% Close Xipmex
xippmex('close');