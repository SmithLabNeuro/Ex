% This example creates stimulation for two different electrodes into one
% command structure.  The first electrode will be stimulated at 30 Hz
% (period = nip clock / Freq = 30kHz/30 = 1000) for one second (repeats =
% desired duration * Freq = 30 * 1 = 30);  The second electrode will have
% all similar stimulation parameters, but will stimulate for two seconds.

% Electrode 1
cmd(1).elec     = 1;
cmd(1).period   = 1000;
cmd(1).repeats  = 30;
cmd(1).action   = 'curcyc';

% Create the first phase (cathodic) for stimulation.  This has a 
% duration of 200 us (6 clock cycles at 30 kHz), an amplitude of 
% 10, and negative polarity.
cmd(1).seq(1) = struct('length', 6, 'ampl', 10, 'pol', 0, ...
    'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', 1);
% Create the inter-phase interval.  This has a duration of 100 us
% (3 clock cycles at 30 kHz).  The amplitude is zero.  The 
% stimulation amp is still used so that the stim markers send by 
% the NIP will properly contain this phase.
cmd(1).seq(2) = struct('length', 3, 'ampl', 0, 'pol', 0, 'fs', 0, ...
    'enable', 0, 'delay', 0, 'ampSelect', 1);
% Create the second, anodic phase.  This has a duration of 200 us 
% (6 cycles at 30 kHz), and amplitude of 10, and positive polarity.
cmd(1).seq(3) = struct('length', 6, 'ampl', 10, 'pol', 1, ...
    'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', 1);


% Electrode 2
cmd(2).elec     = 2;
cmd(2).period   = 1000;
cmd(2).repeats  = 60;
cmd(2).action   = 'curcyc';

% Create the first phase (cathodic) for stimulation.  This has a 
% duration of 200 us (6 clock cycles at 30 kHz), an amplitude of 
% 10, and negative polarity.
cmd(2).seq(1) = struct('length', 6, 'ampl', 10, 'pol', 0, ...
    'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', 1);
% Create the inter-phase interval.  This has a duration of 100 us
% (3 clock cycles at 30 kHz).  The amplitude is zero.  The 
% stimulation amp is still used so that the stim markers send by 
% the NIP will properly contain this phase.
cmd(2).seq(2) = struct('length', 3, 'ampl', 0, 'pol', 0, 'fs', 0, ...
    'enable', 0, 'delay', 0, 'ampSelect', 1);
% Create the second, anodic phase.  This has a duration of 200 us 
% (6 cycles at 30 kHz), and amplitude of 10, and positive polarity.
cmd(2).seq(3) = struct('length', 6, 'ampl', 10, 'pol', 1, ...
    'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', 1);


% Send the stimulation
xippmex('stimseq', cmd);


