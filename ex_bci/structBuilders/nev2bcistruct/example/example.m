%% example script
%% add helper function paths
addpath(genpath('../../nev2bcistruct'))

%% constants
filepath = '';
filename = 'Pe190107_s508a_distanceStabilityBCIcalibMGS_0002';
correctcode = 150;
targoff = 140;
fixoff = 3;
binsize = 0.05;

%% create structure
dat = nev2bcistruct([filename,'.nev']);

%% add in .ns5 data
dat = getNS5Data(dat,filepath,filename,'.ns5');

%% detect dropped start or end trial codes
cleandat = removeBadTrials(dat);

%% get correct trials
correctdat = dat([dat.result]==correctcode);

%% get a task parameter (two methods)
% method 1
angles = driftchoiceextractparam(correctdat,'angle=',2); % set to 2 b/c angle is changed mid-trial in some tasks
% method 2
distances = zeros(length(correctdat),1);
for n = 1:length(distances)
    distances(n) = correctdat(n).params.trial.distance;
end

%% unpack spikes from compressed format to readNEV format
spikes = unpackSpikes(correctdat(1).spikeinfo, correctdat(1).spiketimesdiff,correctdat(1).firstspike,30000);

%% compute spike counts for all trials
% note: many of these outputs are specific to a memory guided saccade task.
% The function will work on other tasks, but these outputs will be empty
[countmatpre,correctdat,] = prepCalibCounts(correctdat,targoff, fixoff,binsize);

%% get psth
% Note: this is not a smoothed psth. If you want it smoothed, then will
% need to apply smoothing kernel to output of the psth2017 function
psth = psth2017(correctdat);

%% get eye trace in pixels and plot it for trial 1
trialnum = 1;
plottraceflag = 1;
eyetrace=getdatEyeTrace(correctdat,trialnum,plottraceflag);