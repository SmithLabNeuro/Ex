% create the overall header values defining electrode, frequency,
% and number of repeats.  This will stimulate on electrode 1 at 30 
% Hz for one second.
cmd = struct('elec', 1, 'period', 1000, 'repeats', 30);
% Create the first phase (cathodic) for stimulation.  This has a 
% duration of 200 us (6 clock cycles at 30 kHz), an amplitude of 
% 10, and negative polarity.
cmd.seq(1) = struct('length', 6, 'ampl', 10, 'pol', 0, ...
    'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', 1);
% Create the inter-phase interval.  This has a duration of 100 us
% (3 clock cycles at 30 kHz).  The amplitude is zero.  The 
% stimulation amp is still used so that the stim markers send by 
% the NIP will properly contain this phase.
cmd.seq(2) = struct('length', 3, 'ampl', 0, 'pol', 0, 'fs', 0, ...
    'enable', 0, 'delay', 0, 'ampSelect', 1);
% Create the second, anodic phase.  This has a duration of 200 us 
% (6 cycles at 30 kHz), and amplitude of 10, and positive polarity.
cmd.seq(3) = struct('length', 6, 'ampl', 10, 'pol', 1, ...
    'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', 1);
% Send the stimulation
xippmex('stimseq', cmd);

