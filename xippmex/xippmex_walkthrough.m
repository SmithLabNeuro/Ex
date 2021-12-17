% This is a Xippmex test script that will run through the following
% operations:
%
% 1) Initialize Xippmex
% 2) Find all Micro+Stim, Micro, Micro HV, and Nano channels (if they exist)
% 3) Configure the data streams on the Front Ends
% 4) Set filters for the Front End LFP and HiRes data streams
% 5) Create and execute a stimulation pattern for a Stim FE
% 6) Collect data
% 7) Plot all collected data included delivered stimulation
% 8) Close Xippmex


%% Initializations 
% Clean the world
close all; fclose('all'); clc; clear all; 

% Initialize xippmex
status = xippmex;
if status ~= 1; error('Xippmex Did Not Initialize');  end

% Find all Stim and Micro/Nano channels and Corresponding FE's
stimChans  = xippmex('elec','stim');
stimFEs = unique(ceil(stimChans/32));

recChans   = [xippmex('elec','micro'), xippmex('elec','nano')];
recFEs = unique(ceil(recChans/32));

% Get NIP clock time right before turning streams on (30 kHz sampling)
timeZero = xippmex('time');

%% Data Stream Activation
% Configure electrode stim FE to have Stim data stream 
% active, all Hi-Res streams to active, and deactivate all other streams
% Note only stim and spike streams are managed individually
tic;
if ~isempty(stimChans)
    xippmex('signal',stimChans,             'stim',   ones(1,length(stimChans)));
    xippmex('signal',stimFEs,               'hi-res', ones(1,length(stimFEs)));   
    xippmex('signal',stimFEs,               'raw',    zeros(1,length(stimFEs)));
    xippmex('signal',stimFEs,               'lfp',    zeros(1,length(stimFEs)));
    xippmex('signal',stimChans,             'spk',    zeros(1,length(stimChans)));
end

% Configure the first 8 recording channels to have Spike streams active
% LFP active on all, and deactivate all other streams
if ~isempty(recChans)
    xippmex('signal',recChans,              'spk',    zeros(1,length(recChans))); 
    xippmex('signal',recChans(1:8),         'spk',    ones(1,8));
    xippmex('signal',recFEs,                'lfp',    ones(1,length(recFEs)));   
    xippmex('signal',recFEs,                'raw',    zeros(1,length(recFEs)));
    xippmex('signal',recFEs,                'hi-res', zeros(1,length(recFEs)));
end
toc;

%% Filter Settings
% Set the Hi-Res notch filter for the first Stim FE to use
% the 60/120/180 comb filter
if ~isempty(stimChans)
    [selection, filters] = xippmex('filter','list',stimFEs(1),'hires notch');
    filterNum = find(strcmp({filters.label}, '60/120/180 Hz'));
    xippmex('filter','set', stimFEs(1),'hires notch', filterNum);
end

% Set the LFP filter for the first recording FE to use the wide-band 1.0-175 Hz filter
if ~isempty(recChans)
    [selection, filters] = xippmex('filter','list',recFEs(1), 'lfp');
    % the below one-liner is simply a trick to which filter label
    % (as cell arrays with strings) have the a low pass of 175 Hz.
    filterNum = find(cellfun(@(x) ~isempty(x), strfind({filters.label}, '175 Hz')));
    xippmex('filter','set', recFEs(1),'lfp', filterNum);
end

%% Pause for 5 seconds to ensure 5 seconds of recording
pause(5); 

%% Stimulation Settings
% Create a stimulation pattern for electrode 5&6 on the first Stim FE with 
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
    
    % Enable stimulation
    xippmex('stim', 'enable', 1);
    
    % Clear Stim Buffer
    xippmex('spike',stimChans(5:6),[1 1]);
    
    % Execute stimulation
    xippmex('stim',stimString);
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
if stimCount(1)
    stims = zeros(3e4*2.5,2);                          % create array of zeros to append stim waveforms to 
    startTime = stimTimestamps{1}(1) - timeZero;     % offset timestamps based on initial NIP clock sample
    
    % loop through each channel, then loop through each stimulation event
    % noted by a stimCount/stimTimestamp. Find when each stim event occured
    % in relation to the first stim event (stimWindow) and append the
    % extracted stimWaveform to the stims array.
    for i = 1:2        
        for j = 1:stimCount(i)
            stimWindow = (1:52) + (stimTimestamps{i}(j) - stimTimestamps{1}(1));         
            stims(stimWindow,i) = stims(stimWindow,i) + ((stimWaveforms{i}(j,:)))';
        end      
    end
    
    time = double(startTime:startTime+size(stims,1)-1)/3e4;  % create array of time for plot usage
    figure(); plot(time,stims(:,1)/1e3, time,stims(:,2)/1e3);
    xlabel('Time (s)'); ylabel('Stim Voltage (mV)'); title('Stimulation Waveforms'); legend('Stim Chan 5','Stim Chan 6'); 
end

%% Plot HiRes Data
if exist('HiResData','var')   
    startTime = floor((HiResTimestamp - timeZero) / 3e4 * 2e3);   % determine when the first HiRes data sample was in relation to the first NIP clock time pulled earlier.
    time = double(startTime:startTime+5*2e3-1)/2e3;                     % create array of time for plot usage
    figure(); plot(time,HiResData(1,:), time,HiResData(2,:), time,HiResData(3,:), time,HiResData(4,:));
    xlabel('Time (s)'); ylabel('Recorded Voltage (uV)'); title('Hi-Res Data'); legend('Stim Ch 9','Stim Ch 10','Stim Ch 11','Stim Ch 12');    
end

%% Plot Spike Data for Channel 3 on first recording FE (if there are any spikes)
if ~isempty(recChans)
    time = (0:51)/30;  % 52 sample waveform at 30kHz sampling
    if spkCount(3)
        figure(); hold on; 
        for i = 1:size(spkWaveforms{3},1); plot(time,spkWaveforms{3}(i,:)); end
        hold off; xlabel('Time (ms)'); ylabel('Recorded Voltage (uV)'); title('Spike Waveforms for Rec Ch 3');
    end
end

%% Plot LFP Data
if ~isempty(recChans)
    startTime = floor((LFPTimestamp - timeZero) / 3e4 * 1e3);  % determine when the first HiRes data sample was in relation to the first NIP clock time pulled earlier.
    time = double(startTime:startTime+5*1e3-1)/1e3;                  % create array of time for plot usage
    figure(); plot(time,LFPData(1,:), time,LFPData(2,:), time,LFPData(3,:), time,LFPData(4,:));
    xlabel('Time (s)'); ylabel('Recorded Voltage (uV)'); title('LFP Data'); legend('Rec Ch 1','Rec Ch 2','Rec Ch 3','Rec Ch 4');    
end


%% Close Xipmex
xippmex('close');