function y = get_stim_word(varargin)
% GET_STIM_WORD - returns stimulation sequence control structure.
% y = get_stim_word(repeats, amp, polarity, enable, 
%                   [fast_settle, delay_length, amp_select])
% returns a structure with the fields 'length', 'amp', 'polarity',
% 'enable', 'fs', 'delay', and 'ampSelect'.
%    LENGTH - number of times for this control word to be repeated (minimum
%       of 1)
%    AMPL - dimensionless stimulation amplitude.  The results stimulation
%       will depend on the configuration of the Micro+Stim frontend.  Contact
%       Ripple to get conversion factors.  Must be 0-127.
%    POL - polarity of stimulation 1 is positive, 0 is negative.  Must be 0
%       or 1
%    ENABLE - Whether is stimulation period transitions to 'on' or 'off'.
%       Must be 0 or 1.
%    FAST_SETTTLE - Whether fast_settle of the neural amp is enabled for
%       this stimulation period.  Must be 0 or 1.
%    DELAY_LENGTH - Length of time in 32nds of a 33.3 us clock cycle (0.04
%       us) to delay this on or off period.  Must be 0 or 1.
%    AMP_SELECT - Selection of neural amp or stimulation amp for this
%       stimulation period.  0 for stimulation amp, 1 for neural amp.
if nargin < 4
    error('usage: get_stim_word(repeats, amp, polarity, enable, [fast_settle, delay_length, amp_select])')
end
% TODO: all these variables need check on bounds and types.
if nargin >= 7
    fast_settle = varargin{5};
    delay_length = varargin{6};
    amp_select = varargin{7};
elseif nargin >= 6
    fast_settle = varargin{5};
    delay_length = varargin{6};
    amp_select = 0;
elseif nargin == 5
    fast_settle = varargin{5};
    delay_length = 0;
    amp_select = 0;
elseif nargin == 4
    fast_settle = 0;
    delay_length = 0;
    amp_select = 0;
end

y = struct('length', varargin{1}, 'ampl', varargin{2}, ...
    'pol', varargin{3}, 'enable', varargin{4}, 'fs', fast_settle, ...
    'delay', delay_length, 'ampSelect', amp_select);
