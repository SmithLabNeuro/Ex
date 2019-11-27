function playTone(freq,length)
% function s = playTone(freq,length)
% freq: frequency in Hz
% length: sound length in seconds

global pahandle;

rate = 48000;
len = round(rate*length);
s = sin((1:len)*2*pi*freq/rate);

s(1:50) = s(1:50).*(1:50)/50;
s(end-49:end) = s(end-49:end).*(50:-1:1)/50;

% Fill the audio playback buffer with the audio data 'wavedata':
PsychPortAudio('FillBuffer', pahandle, s);

% Start audio playback for 'repetitions' repetitions of the sound data,
% start it immediately (0) and wait for the playback to start, return onset
% timestamp.
t1 = PsychPortAudio('Start', pahandle, 1, 0, 1);
