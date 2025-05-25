function alignment_times = bcidemo_alignmentTimesFactory(data, align_to,trials)

% This factory finds the time in the trial in which an event occurs (for
% example by finding a code)
% Inputs:
    % data - data strucute of the type that is created by bcidemo_get_data
    % aligh_to - type of event to align to
    % trials - idices of trials to look for aligment times in
% outputs: alignment_times - vector of the time sin whoch the event happens


global codes
alignment_times = nan(1,length(trials));


switch align_to

    case {'ALIGN','STIM2_ON','STIM1_ON','TARG_ON','STIM_OFF', 'CURSOR_ON','ACQUIRE_TARG'}

        for i = 1:length(trials)
            trial = data.trials(trials(i));           
            alignment_times(i) = align_by_code(trial,codes.(align_to));
        end

    case 'TARG_ONSETS'

        for i = 1:length(trials)
            trial = data.trials(trials(i));           
            alignment_times(i) = align_by_code_set(trial,codes.('TARG_ON'):codes.('TARG10_ON'));
        end

    otherwise
        error('No aligment!')
end
end


function time = align_by_code(trial,code)


inx = find(trial.event_types == code);

assert(length(inx)==1)

time = trial.event_times(inx);

end


function time = align_by_code_set(trial,code_set)


inx = find(ismember(trial.event_types,code_set));

assert(length(inx)==1)

time = trial.event_times(inx);

end