%%% GLOBAL PARAMETERS FILE

%%%% 
debug = 0; % 1 to turn on debug mode, 0 to leave it off (or just comment out these two lines)

%default local folder locations:
if isunix
    rootDir = '~';
else
    error('RUNEX:BadPlatform','Ex only works on unix');
end
params.localDataDir = [rootDir,filesep,'exData'];
params.localExDir = [rootDir,filesep,'Ex_local'];
params.localShowexDir = [rootDir,filesep,'showexLocal'];

%%% settings for randomized rewards:
params.randomizeReward              = false;
params.rewardRandomizationParameter = 'juiceX';
params.rewardDistribution           = 'unif';   %matlab string for the family of distributions from which the reward info will be drawn. See RANDOM for details.
params.rewardDistributionParams     = [1,5]; %Vector containing the parameters for the distribution from which the randomized reward info will be drawn.

%%%% Settings for debugging/demo
params.getEyes = true; % true for using monkey eye movements, false for mouse
params.eyeTrackerAnalog = true; % true for using analog card, false for eyelink via ethernet
params.sendingCodes = true; % sending digital codes
params.rewarding = true; % providing rewards
params.writeFile = true; % write trial data to file
params.bciEnabled = true; % enable BCI computer communication
params.bciCursorEnabled = true; % enable BCI cursor on runex computer
params.statusUpdates = false; % writing a status file of behavioral info
params.soundEnabled = true; % enable sound capabilities
%
params.control2displayIP = '192.168.1.11'; % local IP address of control computer
params.display2controlIP = '192.168.1.10'; % local IP address of display computer
params.control2displaySocket = 4243;
%
params.control2bciIP = '192.168.2.11'; % local IP address of control computer
params.bci2controlIP = '192.168.2.10'; % local IP address of bci computer
params.control2bciSocket = 4244;
%
params.screenDistance = 36; % distance from eye to screen in cm
params.pixPerCM = 27.03; % pixels per centimeter of screen
% fixation window (pixels)
params.fixWinRad = 30; %20  -use a column vector (e.g., [20;20]) for a rectangular (or square) window -first element is half-width and second element is half-height)
% flag to decide whether to recenter the fixation window (in waitForMS)
params.recenterFixWin = false;
% saccade window (pixels) - used when detecting saccade from current eye position
params.sacWinRad = params.fixWinRad/2;
% target window (pixels)
params.targWinRad = 55; %[85;35]; %35
% diode Location (left/top/right/bottom of diode position)
params.diodeLoc = [10 10 40 40];
% if > this time elapses between while loop iterations in waitFor, a warning
% is displayed (to notify us that there are timing issues that may warrant attention)
params.waitForTolerance = 0.01; % in seconds, so 0.01 = 10 ms

%% sound params
params.sampleFreq = 48000; % 48 kHz
params.outBufferSize = floor(params.sampleFreq * 10); % 10 seconds

%% calibration params
params.extent = 250; % spacing of calibration dots in pixels
params.calibX = [-1 0 1] * params.extent;
params.calibY = [1 0 -1] * params.extent;
% juice-related params
params.juiceX = 1; % number of times juice is repeated
params.juiceInterval = 150; % in ms
params.juiceTTLDuration = 30; % in ms, must be >= 1 to give juice
params.juiceChan = 18; % C2, C1/C0 are used for strobe
% other TTL output parameters - C3-C7 are available, bits 19-23
params.digOut1 = 19; % this is C3, typically used for microstim
params.digOut2 = 20; % C4
params.digOut3 = 21; % C5
params.digOut4 = 22; % C6
% used by plotEyes to smooth eye movements - currently just a mean of last
% 'n' data points
params.eyeSmoothing = 5; % 0 for no smoothing, numbers >= 1 are in msec
params.eyeHistoryBufferSize = 200; %in samples (duration depends on how fast "samp" is called)
params.drawSaccades = true;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 DIGITAL CODES LIST (0-255)            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% General instructions for codes
% 0-255: standard trial codes
% 256-511: ascii text
% 1000-32000: custom user code range - use these for exfile-specific codes
% 32768-65535: condition numbers
codes = struct();

% trial boundaries
codes.START_TRIAL = 1;
codes.SHOWEX_TIMINGERROR = 251;
codes.BCI_ABORT = 252;  % use to indicate that BCI computer is not working
codes.ALIGN = 253; % use if you want to align your data to a time point
codes.SHOWEX_ABORT = 254;
codes.END_TRIAL = 255;

% stimulus / trial event codes
codes.FIX_ON = 2 ;
codes.FIX_OFF = 3 ;
codes.FIX_MOVE = 4 ;
codes.REWARD = 5 ;
codes.DIODE_ON = 6 ;
codes.DIODE_OFF = 7 ;
%
codes.NO_STIM = 9;
codes.STIM_ON = 10 ;
codes.STIM1_ON = 11 ;
codes.STIM2_ON = 12 ;
codes.STIM3_ON = 13 ;
codes.STIM4_ON = 14 ;
codes.STIM5_ON = 15 ;
codes.STIM6_ON = 16 ;
codes.STIM7_ON = 17 ;
codes.STIM8_ON = 18 ;
codes.STIM9_ON = 19 ;
codes.STIM10_ON = 20 ;
codes.STIM_OFF = 40 ;
codes.STIM1_OFF = 41 ;
codes.STIM2_OFF = 42 ;
codes.STIM3_OFF = 43 ;
codes.STIM4_OFF = 44 ;
codes.STIM5_OFF = 45 ;
codes.STIM6_OFF = 46 ;
codes.STIM7_OFF = 47 ;
codes.STIM8_OFF = 48 ;
codes.STIM9_OFF = 49 ;
codes.STIM10_OFF = 50 ;
codes.TARG_ON = 70 ;
codes.TARG1_ON = 71 ;
codes.TARG2_ON = 72 ;
codes.TARG3_ON = 73 ;
codes.TARG4_ON = 74 ;
codes.TARG5_ON = 75 ;
codes.TARG6_ON = 76 ;
codes.TARG7_ON = 77 ;
codes.TARG8_ON = 78 ;
codes.TARG9_ON = 79 ;
codes.TARG10_ON = 80 ;
codes.TARG_OFF = 100 ;
codes.TARG1_OFF = 101 ;
codes.TARG2_OFF = 102 ;
codes.TARG3_OFF = 103 ;
codes.TARG4_OFF = 104 ;
codes.TARG5_OFF = 105 ;
codes.TARG6_OFF = 106 ;
codes.TARG7_OFF = 107 ;
codes.TARG8_OFF = 108 ;
codes.TARG9_OFF = 109 ;
codes.TARG10_OFF = 110 ;
codes.CHOICE0 = 120;
codes.CHOICE1 = 121;
codes.CHOICE2 = 122;
codes.CHOICE3 = 123;
codes.CHOICE4 = 124;
codes.CHOICE5 = 125;
codes.CHOICE6 = 126;
codes.CHOICE7 = 127;
codes.CHOICE8 = 128;
codes.CHOICE9 = 129;

% ustim codes
codes.USTIM_ON = 130 ;
codes.USTIM_OFF = 131 ;

% sound codes
codes.SOUND_ON = 132 ;
codes.SOUND_OFF = 133 ;
codes.SOUND_CHANGE = 134;

% behavior codes
codes.FIXATE  = 140 ;	% attained fixation 
codes.SACCADE = 141 ;	% initiated saccade

% trial outcome codes
codes.CORRECT = 150 ;	% Independent of whether reward is given
codes.IGNORED = 151 ;	% Never fixated or started trial
codes.BROKE_FIX = 152 ; % Left fixation before trial complete
codes.WRONG_TARG = 153 ; % Chose wrong target
codes.BROKE_TARG = 154 ; % Left target fixation before required time
codes.MISSED = 155 ;	% for a detection task
codes.FALSEALARM = 156 ;
codes.NO_CHOICE = 157 ;	% saccade to non-target / failure to leave fix window
codes.WITHHOLD = 158 ; % correctly-withheld response
codes.ACQUIRE_TARG = 159 ; % Acquired the target
codes.FALSE_START = 160 ; % left too early
codes.BCI_CORRECT = 161 ; % BCI task (vs. non-BCI behavior) performed correct
codes.BCI_MISSED = 162 ; % BCI task (vs. non-BCI behavior) performed incorrectly 
codes.CORRECT_REJECT = 163 ;
codes.LATE_CHOICE = 164 ;

%%
% retry.CORRECT = 0 ;	% Independent of whether reward is given
% retry.IGNORED = 1 ;	% Never fixated or started trial
% retry.BROKE_FIX = 1 ; % Left fixation before trial complete
% retry.WRONG_TARG = 0 ; % Chose wrong target
% retry.BROKE_TARG = 1 ; % Left target fixation before required time
% retry.MISSED = 0 ;	% for a detection task
% retry.FALSEALARM = 0 ;
% retry.NO_CHOICE = 0 ;	% saccade to non-target / failure to leave fix window
% retry.WITHHOLD = 0 ; %correctly-withheld response
% retry.SACCADE = 0;
% touch bar / lever / button press codes would go here

% OLD CODES
% codes for sending over digital port
%codes = struct();
%codes.START_TRIAL = 1;
%codes.END_TRIAL = 2;
%codes.FIX_ON = 5;
%codes.FIX_OFF = 6;
%codes.STIM_ON = 7;
%codes.STIM_OFF = 8;
%codes.FIX_MOVE = 9;
%codes.FIX_CAUGHT = 10;
%codes.CONDITION = 15;
%codes.JUICE = 18;
%codes.REWARD = 19;
%codes.START_HISTOGRAM = 50;
%codes.ALIGN_HISTOGRAM = 51;
%codes.STOP_HISTOGRAM = 52;
% codes.NOFIXATION
% codes.BROKEFIXATION (during the stimulus)
% codes.DIDNTGETTOTARGET
% codes.LEFTTARGETEARLY
% codes.WRONGTARGET
% codes.CORRECTTARGET
% general ABORT code?

%% assign the value of 'debug' in runex.m:
assignin('caller','debug',debug);

