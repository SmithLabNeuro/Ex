function data = bcidemo_getData(dat)

% This function organizes the data for the other BCI functions
% The input is a dat structure, the output has a similar format but is
% compatible with the BCI fucntions in the utils


for t = 1:length(dat)

    data.trials(t).spike_times = get_spike_times([dat(t).firstspike], ...
        [dat(t).spiketimesdiff], [dat(t).spikeinfo]);
    data.trials(t).spike_neuro_ids = get_neuron_ids_for_spikes([dat(t).spikeinfo]);
    data.trials(t).event_types = [dat(t).trialcodes(:,2)];
    data.trials(t).event_times = [dat(t).trialcodes(:,3)];
    data.trials(t).outcome = trial_outcome(dat(t));
    data.trials(t).trial_onset =  find_trial_onset(dat(t));
    data.trials(t).trial_offset =  find_trial_offset(dat(t));
    data.trials(t).target = find_target_direction(dat(t));
end

data.info.all_ids = unique([data.trials.spike_neuro_ids]);
end

function spike_times = get_spike_times(first_spike, spike_diff, channels)

SAMPLING_RATE = 30000;

% Get spike times from raw spike data.

% Parameters:
% - first_spike (array): Array containing first spike times.
% - spike_diff (array): Array containing spike time differences.
% - channels (array): Array containing channel information.

% Returns:
% - spike_times (array): Array containing spike times.
spike_diff = spike_diff';
spike_times = [first_spike, first_spike+cumsum(double(spike_diff))];
spike_times = spike_times / SAMPLING_RATE;
end

function neuron_ids = get_neuron_ids_for_spikes(channels)
% Get neuron IDs corresponding to spikes.

% Parameters:
% - channels (array): Array containing channel information.

% Returns:
% - neuron_ids (array): Array containing neuron IDs.

neuron_ids = channels(:,1)';

end


function outcome = trial_outcome(trial)

global codes
trial_codes = [trial.trialcodes(:,2)];

POSSIBLE_OUTCOMES = setdiff([codes.BACKGROUND_PROCESS_TRIAL, codes.CORRECT:codes.BROKE_TASK],[codes.ACQUIRE_TARG codes.WITHHOLD]);

inx_outcome = find(ismember(trial_codes,POSSIBLE_OUTCOMES));

assert(length(inx_outcome)==1)

outcome = trial_codes(inx_outcome);

end

function completed = trial_completed(trial)

global codes

COMPLETE_OUTCOMES = [codes.CORRECT, codes.CORRECT_REJECT, ...
    codes.FALSEALARM, codes.MISSED,codes.WRONG_TARG,codes.NO_CHOICE];

outcome = trial_outcome(trial);

completed = ismember(outcome,COMPLETE_OUTCOMES);


end

function angle = target_orientation_change(trial)

if trial.params.trial.isCatch
    angle = 0;

else
    angle = abs(trial.params.trial.targAmp);
end

end


function onset = find_trial_onset(trial)

global codes

inx = find(trial.trialcodes(:,2)== codes.START_TRIAL);
assert(length(inx)==1)
onset = trial.trialcodes(inx,3);
end


function onset = find_trial_offset(trial)

global codes

inx = find(trial.trialcodes(:,2)== codes.END_TRIAL);
assert(length(inx)==1)
onset = trial.trialcodes(inx,3);
end


function direction = find_target_direction(trial)
direction = trial.params.trial.targetAngle;
end