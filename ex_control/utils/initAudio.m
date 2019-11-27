function initAudio()

global pahandle;

% Perform basic initialization of the sound driver:
InitializePsychSound;

PsychPortAudio('Close');

% Open the default audio device [], with default mode [] (==Only playback),
% and a required latencyclass of zero 0 == no low-latency mode, as well as
% a frequency of freq and nrchannels sound channels.
% This returns a handle to the audio device:
pahandle = PsychPortAudio('Open', [], [], 1, 48000, 1);