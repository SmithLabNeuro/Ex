function [stim_cmd] = xippmexStimCmd(stim_chan,pulse_width,stim_freq,stim_dur,stim_amp)
%function [stim_cmd] = xippmexStimCmd(stim_chan,pulse_width,stim_freq,stim_dur,stim_amp)
%
% This function formats the stim command that gets sent to Ripple via
% xippmex

% eye movement = 50 uA, 70ms, 350 Hz, 250 uS, negative first?

% 11Jun2015 by ACS: added some support for multiple channels. Right now,
% variable amplitudes can be specificed, but other parameters (pulse width,
% frequency, etc.) are yoked across channels. -ACS 11Jun2015

% CHECKING FOR REASONABLE LIMITS ON USTIM
if pulse_width ~= 250
    error('Why mess with this? Just use 250 for pulse_width');
end
if stim_freq ~= 350
    error('Why mess with this? Just use 350 for stim_freq');
end
if stim_dur <= 0 || stim_dur > 500
    error('Are you trying to fry things? Keep stim_dur between 0 and 500');
end
if any(stim_amp < 0 | stim_amp > 250)
    error('Are you trying to fry things? Keep stim_amp between 0 and 250');
end
if ~any(stim_amp)
    error('At least one channel must have a non-zero stim_amp');
end

if numel(stim_amp) > 1
    if numel(stim_chan) ~= numel(stim_amp)
        error('Must specify either 1 amp for all channels or 1 amp for each channel');
    end
end

clock_cycle = 1/30 * 1000; % clock ticks in microseconds
delay_length = clock_cycle / 32; % length of a single unit of delay in microseconds

% basic FEF eye movement-inducing parameters
%stim_chan = 1; % channel number to stimulate on
%pulse_width = 250; % uS of each phase (cathodic or anodic)
%stim_freq = 350; % Hz of stimulation pulses
%stim_dur = 100; % ms of stimulation
%stim_amp = 100; % uA of stimulation

% for micro+stim it's hard set to 7.5, for micro2/nano2+stim it's
% adjustable
%1 = 1 uA/step, 2 = 2 uA/step, 3 = 5 uA/step, 4 = 10 uA/step, 5 = 20 uA/step
%step_size = 7.5; % uA of each step built into the micro+stim

% loop and set 5uA/step for all channels
disableTimeout = 1;
for I=1:numel(stim_chan)
    disableStart = tic;
    while xippmex('stim','enable') %use a while loop to confirm that this has actually been disabled -acs26apr2016
        assert(toc(disableStart)<disableTimeout,'xippmexStimCmd:disableTimeout','Timeout disabling stimulation to set step size');
        xippmex('stim','enable',0);
    end;
    xippmex('stim','res',stim_chan(I),3);
end
step_size = 5.0; % uA of each step built into the micro+stim

% convert our stim duration/freq to what the micro+stim wants
stim_clockticks = (stim_dur * 1000) / clock_cycle;
stim_repeats = floor((stim_dur / 1000) * stim_freq);
stim_period = floor(stim_clockticks / stim_repeats);
full_pulses = floor(pulse_width / clock_cycle);

% figure out step size for the micro+stim
nstim_steps = floor(stim_amp/step_size);
%nstim_steps(nstim_steps<1) = 1; %changed from a conditional to support vector arguments -ACS 11Jun2015
actual_stim_amp = nstim_steps * step_size; % actual stimulation amplitude due to discrete steps

% get rid of stim_channels where the amp is zero - needs to be updated if
% we allow variable pulse widths, etc
idx = find(stim_amp==0);
if ~isempty(idx)
    actual_stim_amp(idx) = [];
    stim_amp(idx) = [];
    stim_chan(idx) = [];
    nstim_steps(idx) = [];
end

disp(['XIPPMEX: Requested ',num2str(stim_amp),' uA and actually delivered ',num2str(actual_stim_amp),' uA on channels ', num2str(stim_chan)]); %might want to add channel numbers to this warning -ACS 11Jun2015

% setup the overall stim command structure
stim_cmd = struct('elec', num2cell(stim_chan), 'period', num2cell(stim_period), 'repeats', num2cell(stim_repeats)); %put arguments in num2cell calls for multiple channel support -ACS 11Jun2015
if isscalar(nstim_steps)&&~isscalar(stim_cmd),
    nstim_steps = nstim_steps.*ones(size(stim_cmd)); %set all amplitudes equal if multiple channels are requested but only one amplitude is supplied. -ACS 11Jun2015
end;
for thisChan = 1:numel(stim_cmd),
    %%% LOOKS LIKE THIS WAS POSITIVE FIRST, WE WANT NEGATIVE FIRST SO I CHANGED
    %%% THE POLARITY ORDER SO POL OF 1 WAS FIRST AND 0 SECOND (DIFFERENT FROM
    %%% THE DEFAULT IN THE RIPPLE CODE)
    %%% 11Jun2015- update on comment above: this seems to have been due to
    %%% the "invert" setting being on on the oscilloscope. I switched the
    %%% polarity order back -ACS 11Jun2015
    
    % setup the cathodic phase - 1 clock cycle for the first 33.3 uS, then fill
    % in the rest of the pulse to get exactly the length you want (with near-uS
    % precision)
    stim_cmd(thisChan).seq(1) = struct('length', full_pulses, 'ampl', nstim_steps(thisChan), 'pol', 0,'fs', 1, 'enable', 1, 'delay', 0, 'ampSelect', 1);
    cath_remaining = pulse_width - (full_pulses * clock_cycle);
    cath_delay = floor(cath_remaining / delay_length);
    stim_cmd(thisChan).seq(2) = struct('length', 1, 'ampl', nstim_steps(thisChan), 'pol', 0,'fs', 1, 'enable', 0, 'delay', cath_delay, 'ampSelect', 1);
    
    % could setup an interphase interval here, but we don't usually use one
    % MATT - optional here, this code as is will wait for the end of the clock
    % cycle to start the anodic phase. So it will really be 248.1 uS instead of
    % 250 uS (for example).
    
    % setup the anodic phase
    stim_cmd(thisChan).seq(3) = struct('length', full_pulses, 'ampl', nstim_steps(thisChan), 'pol', 1,'fs', 1, 'enable', 1, 'delay', 0, 'ampSelect', 1);
    an_remaining = pulse_width - (full_pulses * clock_cycle);
    an_delay = floor(an_remaining / delay_length);
    stim_cmd(thisChan).seq(4) = struct('length', 1, 'ampl', nstim_steps(thisChan), 'pol', 1,'fs', 1, 'enable', 0, 'delay', an_delay, 'ampSelect', 1);
    
end;

end


%% SCRAP COMMENTS FROM THE DEVELOPMENT PHASE:
% simple command
%stim_cmd.seq(1) = struct('length', 8, 'ampl', nstim_steps, 'pol', 0,'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', 1);
%stim_cmd.seq(2) = struct('length', 8, 'ampl', nstim_steps, 'pol', 1,'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', 1);
%return

% stim_cmd = struct('elec', stim_chan, 'period', 1, 'repeats', 1);
% stim_cmd.seq(1) = struct('length', 6, 'ampl', nstim_steps, 'pol', 0,'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', 1);
% stim_cmd.seq(2) = struct('length', 3, 'ampl', 0, 'pol', 0,'fs', 0, 'enable', 0, 'delay', 0, 'ampSelect', 1);
% stim_cmd.seq(3) = struct('length', 6, 'ampl', nstim_steps, 'pol', 1,'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', 1);
% return

%microstim parameters
% create the overall header values defining electrode, frequency,
% and number of repeats. This will stimulate on electrode 1 at 200
% Hz for 50 ms.
%cmd = struct('elec', 60, 'period', 150, 'repeats', 10);
% Create the first phase (cathodic) for stimulation. This has a
% duration of 100 us (3 clock cycles at 30 kHz), an amplitude of
% 5, and negative polarity.
%cmd.seq(1) = struct('length', 3, 'ampl', 5, 'pol', 0, ...
%'fs', 1, 'enable', 1, 'delay', 0, 'ampSelect', 1);
% Create the second, anodic phase. This has a duration of 100 us
% (3 cycles at 30 kHz), and amplitude of 5, and positive polarity.
%cmd.seq(2) = struct('length', 3, 'ampl', 5, 'pol', 1, ...
%'fs', 1, 'enable', 1, 'delay', 0, 'ampSelect', 1);


% % %second channel
%cmd1 = struct('elec', 21, 'period', 150, 'repeats', 10);
% Create the first phase (cathodic) for stimulation. This has a
% duration of 100 us (3 clock cycles at 30 kHz), an amplitude of
% 5, and negative polarity.
%cmd1.seq(1) = struct('length', 3, 'ampl', 5, 'pol', 0, ...
%'fs', 1, 'enable', 1, 'delay', 0, 'ampSelect', 1);
% Create the second, anodic phase. This has a duration of 100 us
% (3 cycles at 30 kHz), and amplitude of 5, and positive polarity.
%cmd1.seq(2) = struct('length', 3, 'ampl', 5, 'pol', 1, ...
%'fs', 1, 'enable', 1, 'delay', 0, 'ampSelect', 1);

%
%
% % % %third channel
% cmd2 = struct('elec', 82, 'period', 150, 'repeats', 10);
% % Create the first phase (cathodic) for stimulation. This has a
% % duration of 100 us (3 clock cycles at 30 kHz), an amplitude of
% % 5, and negative polarity.
% cmd2.seq(1) = struct('length', 3, 'ampl', 5, 'pol', 0, ...
% 'fs', 1, 'enable', 1, 'delay', 0, 'ampSelect', 1);
% % Create the second, anodic phase. This has a duration of 100 us
% % (3 cycles at 30 kHz), and amplitude of 5, and positive polarity.
% cmd2.seq(2) = struct('length', 3, 'ampl', 5, 'pol', 1, ...
% 'fs', 1, 'enable', 1, 'delay', 0, 'ampSelect', 1);
%
% %stim_cmd=[cmd];
% %stim_cmd=[cmd, cmd1];
% stim_cmd=[cmd, cmd1, cmd2];
%

% Send the stimulation
%xippmex('stimseq', stim_cmd);