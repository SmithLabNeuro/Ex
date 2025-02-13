function raster = bcidemo_calculateRasters(data, align_to,time_before,time_after,trials,id)

% claculate raster for a specific unit. 
% data: Structured data from bcidemo_getData.
% align_to: Specifies what event or time point to align the spikes to.
% time_before: Time window before the alignment event.
% time_after: Time window after the alignment event.
% trials: Indices of trials to process.
% id: ID of the neuron whose spikes are being analyzed.
% Output:
% raster: A cell array containing raster  for each specified trial.

alignment_times = bcidemo_alignmentTimesFactory(data,align_to,trials);
raster = cell(1,length(trials));

for i=1:length(trials)
    spike_times = data.trials(trials(i)).spike_times;
    spike_times = spike_times - alignment_times(i);
    spike_ids = data.trials(trials(i)).spike_neuro_ids;
    rel_spikes = (spike_times > time_before) & ((spike_times < time_after) ...
       & (spike_ids == id));
    raster{i} = spike_times(rel_spikes);
end


