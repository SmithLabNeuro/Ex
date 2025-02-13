function ff = bcidemo_getChanFF(data,chan)

% Inputs:
% data: Structured data containing information about trials, 
% including trial onset, trial offset, and spike times.
% chan: Channel number for which the firing rate needs to be calculated.
% Output:
% ff: The ff of the specified channel

rate_per_trial = nan(1,length(data.trials));

for t = 1:length(data.trials)   
    duration = (data.trials(t).trial_offset - data.trials(t).trial_onset);
    spikes = bcidemo_getTrialSpikesByChan(data.trials(t),chan);
    rate_per_trial(t) = length(spikes)/duration;
end

ff = var(rate_per_trial)/mean(rate_per_trial); 