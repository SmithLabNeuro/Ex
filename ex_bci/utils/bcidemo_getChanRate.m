function rate = bcidemo_getChanRate(data,chan)

% Inputs:
% data: Structured data containing information about trials, 
% including trial onset, trial offset, and spike times.
% chan: Channel number for which the firing rate needs to be calculated.
% Output:
% rate: The firing rate of the specified channel, expressed in Hz 
% (spikes per second).

duration = 0;
num_spikes = 0;
for t = 1:length(data.trials)   
    duration = duration + (data.trials(t).trial_offset - data.trials(t).trial_onset);
    spikes = bcidemo_getTrialSpikesByChan(data.trials(t),chan);
    num_spikes = num_spikes + length(spikes);
end

rate = num_spikes/duration; % Hz