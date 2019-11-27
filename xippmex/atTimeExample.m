% This example uses the 'at-time' command for stimulation to make a 10
% pulse, 50 Hz stimulation train, with 100 us biphasic stim, execute
% exactly 200 ms from the requested, current NIP processor time, on the 
% first found stimulation channel. This example can be used to create 
% stimulation events that are timed based on spike event timestamps.

% Clean the world
close all; fclose('all'); clc; clear all; 

% Initialize xippmex
status = xippmex;
if status ~= 1; error('Xippmex Did Not Initialize');  end

% Enable stimulation on the NIP.
% NOTE: If stimulation is not enabled, xippmex will enable and disable 
% stimulation for each stim sequence it receives. This will have the 
% undesirable side effect of killing queued stimulation sequences. Always 
% enable stimulation on the NIP before queueing multiple stimulation
% sequences.

xippmex('stim', 'enable', 0); pause(0.5)
xippmex('stim', 'enable', 1); pause(0.5)

% when 'at-time' is used, the cmd.period field changes meaning and 
% represents the lower 16 bits of the NIP clock when the stim should be 
% executed. Only a single pass through the stim sequence is executed, i.e.,
% repeats is not used, but is still required for the structure.
stimChans = xippmex('elec','stim');

cmd = struct('elec',    stimChans(1), ...
             'period',  1, ...
             'repeats', 1, ...
             'action',  'at-time');
         
% Valid 'action' options are: ['immed', 'curcyc', 'allcyc', 'trigger', 'at-time']

% Set stimulation waveform to be 10 pulses at 50 Hz where the waveform has
% phase durations of 100 us (3 clock cycles), biphasic squarewaves, with 
% 33 us (1 clock cycle) inter-phase interval.
pw_us        = 100;
ipi_us       = 100/3;
numPulses    = 10;
stim_freq_Hz = 50;
nipClock_us  = 1e6/3e4;
pw_cycs      = floor(pw_us / nipClock_us);
ipi_cycs     = floor(ipi_us / nipClock_us);
period_cycs  = floor(1e6 / stim_freq_Hz / nipClock_us) - 2*pw_cycs - ipi_cycs;

seqIdx = 1;
for i = 1:numPulses
    
    % Set cathodic phase of stim
    cmd.seq(seqIdx) = struct('length', pw_cycs, 'ampl', 0, 'pol', 0, ...
    'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', 1);
    seqIdx = seqIdx+1;

    % Set interphase interval, i.e., time between cathodic an anodic 
    % phases of bipolar waveform. In this example, it is one clock cycle.
    cmd.seq(seqIdx) = struct('length', ipi_cycs, 'ampl', 0, 'pol', 0, ...
        'fs', 0, 'enable', 0, 'delay', 0, 'ampSelect', 1);
    seqIdx = seqIdx+1;

    % Set anodic phase of stim
    cmd.seq(seqIdx) = struct('length', pw_cycs, 'ampl', 0, 'pol', 1, ...
        'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', 1);
    seqIdx = seqIdx+1;
    
    % Set time between pulses, where stim level is 0. Do not add command
    % sequence here for the last pulse. Ensure all values except length are
    % zero, otherwise NIP/Trellis report stim events during OFF part of
    % period
    %if i ~= numPulses
        cmd.seq(seqIdx) = struct('length', period_cycs, 'ampl', 0, 'pol', 0, ...
            'fs', 0, 'enable', 0, 'delay', 0, 'ampSelect', 0);
        seqIdx = seqIdx+1;
    %end
end
    
% Generate the mask that will be used to set the execution time of each stim seqence
atTimeMask = uint32(hex2dec('FFFF'));

% Get the current NIP time
curNipTime = xippmex('time');

% Set the stimulation execution time (cmd.period) to be 200 ms from the current time
cmd.period = bitand( uint32(curNipTime + floor(200e3 / nipClock_us) ), atTimeMask);

% Send stimulation command
xippmex('stimseq', cmd);

% close xippmex
%xippmex('close');