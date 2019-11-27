function [ y ] = get_stim_cmd(elec, repeats, period)
% GET_STIM_CMD returns stim control structure with proper format
% y = get_stim_cmd(elec, repeats, period) returns a structure with the
% fields 'elec', 'repeats', and 'period'
%   ELEC - one indexed electrode id.
%   REPEATS - number of time to repeat stim pattern.
%   PERIOD - period of stimulation 33.3 us units (30 kHz clock cycles).

% TODO: need to check the bounds of the given parameters here.
y = struct('elec', elec, 'repeats', repeats, ...
        'period', period, 'seq', []);
end

