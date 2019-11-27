% xippmex('stimseq') example producing a symetric biphasic pulse with
% phases of 50 us and an inter-phase interval of 100 us.  This example
% makes use of the stimulation control word's 'delay' field in order to
% produce phases with a duration that is smaller than 33.3 us.

% single clock cycle in micro-seconds
clock_cycle = 1 / 30 * 1000;
% calculate the length of a single unit of delay in micro-seconds
delay_length = 1 / 30 * 1000 / 32;
% Setup the basic header elements.  This will produce stimulation at 30 Hz
% for one second
cmd = struct('elec', 1, 'period', 1000, 'repeats', 30);
% first command word.  Produce 1 clock_cycle worth of stimulation, this
% 33.3 us of the first 50 us phase.
cmd.seq(1) = struct('length', 1, 'ampl', 10, 'pol', 0, 'fs', 0, ...
    'enable', 1, 'delay', 0, 'ampSelect', 1);
% calculate the remaining duration of the cathodic phase.  Because 
% this is less than one clock cycle period the delay field will need
% to be used
cath_remaining = 50.0 - clock_cycle;
% the actual delay parameter must be an integer between 0 and 31
cath_delay = floor(cath_remaining / delay_length);
% This word will produce the full cathodic phase and a bit of the
% inter-phase interval.  Here, the 'enable' field set to zero, sets 
% the pulse to start high and transition to off.
cmd.seq(2) = struct('length', 1, 'ampl', 10, 'pol', 0, 'fs', 0, ...
    'enable', 0, 'delay', cath_delay, 'ampSelect', 1);
% The inter-phase interval is 100 us, so add two more full clock 
% cycles worth of stimulation for 66.6 us more.
cmd.seq(3) = struct('length', 2, 'ampl', 0, 'pol', 0, 'fs', 0, ...
    'enable', 0, 'delay', 0, 'ampSelect', 1);
% calculate how much more of the inter-phase interval is remaining 
% being careful to account for the quantization of the delay at the 
% end of the cathodic phase
ipi_remaining = ...
    100.0 - (clock_cycle - cath_delay * delay_length) - ...
    2 * clock_cycle;
ipi_delay = floor(ipi_remaining / delay_length);
% add this to the sequence list.  Here the enable field sets the
% pulse to start off and transition to on.
cmd.seq(4) = struct('length', 1, 'ampl', 10, 'pol', 1, 'fs', 0, ...
     'enable', 1, 'delay', ipi_delay, 'ampSelect', 1);
% At this point there is not a full clock cycle's worth of stim left,
% so the pulse will be finished with one use of the delay field
 anod_remaining = ...
    50.0 - (clock_cycle - ipi_delay * delay_length);
anod_delay = floor(anod_remaining / delay_length);
% add the last word
cmd.seq(5) = struct('length', 1, 'ampl', 10, 'pol', 1, 'fs', 0, ...
     'enable', 0, 'delay', anod_delay, 'ampSelect', 1);
% fire off stimulation
xippmex('stimseq', cmd);

