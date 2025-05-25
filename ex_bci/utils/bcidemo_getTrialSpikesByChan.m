function spikes = bcidemo_getTrialSpikesByChan(trial,chan)

% This function extracts the spike times from a specific channel
% within a given trial
% 
% Parameters:
% 
% trial: A structure containing the spike data for a single trial. 
% It must include the fields spike_neuro_ids (neuron/unit identifiers) 
% and spike_times (the corresponding spike times).
% chan: channel number
% 
% Returns:
% 
% spikes: A vector of spike times (in seconds or the relevant time unit)
% for the channel identified by id in the specified trial. If no spikes 
% are found for the given ID, it will return an empty vector.

inx = find(trial.spike_neuro_ids==chan);
spikes = trial.spike_times(inx);

end
