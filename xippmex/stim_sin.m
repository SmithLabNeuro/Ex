 function varargout = stim_sin(elec, freq, wf_amp, varargin)
% STIM_SIN produce sinusiodal stimulation.
% Stimulation output will have a sinusoidal resolution of 33.3 us (30 kHz
% clock cycles).  The output range will be discritized in 127 steps in 
% both the positive an negative directions (254 steps total).
% STIM_SIN(ELEC, FREQ, AMP[, REPEATS])
%    ELEC - desired electrode
%    FREQ - frequency of sin wave in Hz 
%    AMPL - sin amplitude in stim step sizes (actual output depends on
%       headstage).
%    REPEATS - number of repeats for entire stimulation pattern

% TODO: check that elec is a valid stim electrode

% state clock freq in Hz
CLOCK_RATE = 30000; 
% number of allowed stim control words allowed in the NIP.
MAX_BUFFER_SIZE = 240;
MAX_REPEATS = hex2dec('fff');
% We'll repeat the train for the maximal amount of time.  The
% overall length however will depend on the input frequency.
if nargin > 3
    repeats = varargin{1};
    if repeats > MAX_REPEATS
        warning('repeats must be less than %d\n', MAX_REPEATS);
        repeats = MAX_REPEATS;
    end
else
    repeats = MAX_REPEATS;
end
% clock counts per period.  Round this to an integer number of clock 
% counts.  This is the input to the NIP command word 'period' must be less 
% than 0xffff.
cc_per_period = floor(CLOCK_RATE / freq);
if cc_per_period > hex2dec('ffff')
    error('desired period is too slow, must be less than %.4f Hz', ...
        hex2dec('ffff') / CLOCK_RATE);
end
if rem(CLOCK_RATE, freq)
    warning('rounding frequency to state clock: %d Hz', ...
        CLOCK_RATE / cc_per_period);
end
if floor(freq) ~= freq
    warning('untest at frequency resolutions less than 1 Hz');
end
if wf_amp < 0 || wf_amp > 120
    error('ensure amplitude between 0 and 120');
end
% create a sine wave in 30kHz bin sizes
xx = linspace(0, 2*pi, cc_per_period);
wf_amp = floor(wf_amp);
yy = round(wf_amp*sin(xx));

cmd = get_stim_cmd(elec, repeats, cc_per_period);
% begin array of structures for stim sequences
seq = [];
% start saved y's at an impossibly high value.  Honestly, this could be
% anything highier than 128
prev_y = 1e9;
% may as well initalize bin_count here
bin_count = 0;
% create a sequence of stim words that combines repeated values in yy using
% the 'length' parameter.  This structure should have the equivalent
% information as the yy array and will be shorter.
for j=1:length(yy)
    % on start, get the first value
    if prev_y == 1e9
        prev_y = yy(j);
        bin_count = 1;
    end
    % if we have a new word or are at the end write the word to the array
    if yy(j) ~= prev_y || j == length(yy)
        ampl = abs(prev_y);
        pol = 1;
        if prev_y < 0
            pol = 0;
        end
        seq = [seq, get_stim_word(bin_count, ampl, pol, 1)];
        % values can be calculated after the fact with an arrayfun, if this
        % ends up being slow.
        
        % ready to go with the search for the next stim word sequence
        prev_y = yy(j);
        bin_count = 1;
    else 
        % since we have another match, iterate the bin_count
        bin_count = bin_count + 1;
    end    
end
% Initially was conserned that the end points wouldn't always line up.  In
% practice this does not seem to be an issue.
% seq(1).ampl = seq(end).ampl;

% if the current words is bigger than the maximum buffer size we'll need
% to trim the size of the words that we're sending
seq_indices = 1:length(seq);
% create a sequence of indices of equal length to remove from the sequece.
% We will need to add the length to the nieghboring bins
indices = round(linspace(1, length(seq), length(seq) - MAX_BUFFER_SIZE));
% throw warning.  Remove for demos?
if ~isempty(indices)
    warning('trimming words');
end
% FIXME: how to handle repeated indices?  Though in practice this doesn't
% seem to be an issue.  May require closer study.
% We start moving backwards so that moved 'lengths' are not lost
for jj=fliplr(1:length(indices))
    removed_index = indices(jj);
    if removed_index == 1
        % skip the first section.  This seems to work fine
        continue;
    end
    % add the length of the current index to that of the previous
    seq(removed_index - 1).length = seq(removed_index-1).length + ...
        seq(removed_index).length;
    % alternate which 'ampl' we keep each time.  I was conserned about
    % systematic effects of pushing to the left sine wave values.  This may
    % be silly in practice.
    if mod(jj, 2) == 0
        seq(removed_index - 1).ampl = seq(removed_index - 1).ampl;
    end
    % store the removed indices
    seq_indices(removed_index) = 0;
end
% remove unwanted sequence words
seq_indices = seq_indices ~= 0;
seq = seq(seq_indices);
% set the seq array in the command word
cmd.seq = seq;
% if no arguments, just do stim
if nargout == 0
    xippmex('stimseq', cmd);
end
% if output is requested don't do stim, but return the command words
if nargout >= 1
    varargout{1} = cmd;
end
if nargout >=2
    varargout{2} = y;
end
