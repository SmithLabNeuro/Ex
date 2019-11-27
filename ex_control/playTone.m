function playTone(freq,duration,speakerind,localizeflag,azimuth,elevation)
%playTone(freq,duration)
%
% Produces a single tone with low latency (~25ms on Macbook Air) using the PsychToolBox
%
% Inputs
% freq: frequency of tone in Hertz
% duration: duration of tone in milliseconds (duration*samplefreq must be
% smaller than buffersize used in initsounds)
% speakerind: takes on value 1 or 2, 1 indicating right and 2 indicating left.
%
% History:
% Written 4/27/2016 Ryan Williamson
% Update 1/23/2017 - added ability to select speaker
% Update 1/24/2017 - added 3D audio

global params audioHandle;

duration = duration / 1000; % convert to seconds for code below

t=linspace(0,duration,duration*params.sampleFreq);

if ~exist('speakerind','var')
    y = (sin(freq*2*pi*t));
    dataout = [y;y];
elseif speakerind == 1 || speakerind == 2
    y = (sin(freq*2*pi*t));
    dataout = zeros(2,length(y));
    dataout(speakerind,:) = y;
elseif speakerind == 3
    if params.sampleFreq ~= 44100;
        error('SampleFreq must equal 44100')
    end
    if ~exist('localizeflag','var')||localizeflag==0
        y = (sin(freq*2*pi*t));
        dataout = [y;y];
    else
        ypre = (sin(freq*2*pi*t));
        load hrir_final.mat
        [pulseleft, ~, ~] = getNearestUCDpulse(azimuth,elevation,hrir_l);
        yl = filter(pulseleft,1,ypre);
        [pulseright, ~, ~] = getNearestUCDpulse(azimuth,elevation,hrir_r);
        yr = filter(pulseright,1,ypre);
        dataout = [yl;yr];
    end
else
    error('Invalid Speaker Index')
end
PsychPortAudio('FillBuffer', audioHandle, dataout);
PsychPortAudio('Start', audioHandle, 1, 0, 1);

end