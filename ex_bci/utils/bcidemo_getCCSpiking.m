function cc = bcidemo_getCCSpiking (data, focus_chan, all_chans, cc_time)

% This function computes the coincident rate between the spikes of a
% focus channel and multiple other channels over a series of trials, 
% based on a specified time window (cc_time).
% 
% Parameters:
% 
% data: Structured data from bcidemo_get_data.
% focus_chan: The index or ID of the focus channel (the channel of interest)
% from which spikes will be compared to other channels.
% all_chans: An array of indices or IDs for all channels that should be 
% compared to the focus channel.
% cc_time: A scalar specifying the maximum time difference (in the same 
% unit as spike times) for two spikes to be considered as "coincident" (ms).

% Returns:
% 
% cc: A row vector of length length(all_chans), where each element
% represents the fraction of trials in which a spike in the focus channel 
% had a corresponding spike in one of the all_chans within the time window 
% cc_time.

cc = zeros(1, length(all_chans));
num_spikes = 0;

for t = 1:length(data.trials)

    focus_spikes = bcidemo_getTrialSpikesByChan(data.trials(t),focus_chan);

    for i = 1:length(all_chans)
        
        other_spikes = bcidemo_getTrialSpikesByChan(data.trials(t),all_chans(i));
        
        diff_mat = abs(focus_spikes' - other_spikes);
        
        spike_has_cc = any(diff_mat < cc_time,2);
        
        cc(i) = cc(i) + sum(spike_has_cc);

    end

    num_spikes = num_spikes + length(focus_spikes);
end

cc = cc / num_spikes; % so that cc will be the fraction of spikes


