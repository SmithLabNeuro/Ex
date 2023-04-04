function [decoderFileLocationAndName] = calibrateInternalStateAxisDecoder(socketsControlComm, nevFilebase, nevFilesForTrain, trainParams, subject)

global params codes

% NASNet variables
gamma = trainParams.gamma;
nasNetName = trainParams.nasNetwork;
netFolder = params.nasNetFolderDataComputer;
% Specify channels that will be used for BCI
channelNumbersUse = trainParams.rippleChannelNumbersInBci;

% don't double read the base if it's also a train file
if any(strcmp(nevFilesForTrain, nevFilebase))
    includeBaseForTrain = true;
else
    includeBaseForTrain = false;
end
nevFilesForTrain(strcmp(nevFilesForTrain, nevFilebase)) = [];

%% Read in NEV files that will be used for training our decoder
[nevBase,waves] = readNEV(nevFilebase);
nev = nevBase;
kr = [];

% Needed to globally align times for multiple NEV files
firstTimePtShift = 0;
fnStartTimes = firstTimePtShift;
% Read in Nev Files and NS5 files that will be used for training decoder
for nevFlInd = 1:length(nevFilesForTrain)
    currNevFileName = nevFilesForTrain{nevFlInd};
    [nevNx, wavesNx] = readNEV(currNevFileName);
    kr = [kr (nev(end, 3) + 1)*30000];
    % Update first Time point based on what nevNx last timestamp is; set to
    % set to 1 second after current nev file's final time stamp
    firstTimePtShift = firstTimePtShift + nevNx(end, 3) + 1;
    fnStartTimes(end+1 ) = firstTimePtShift;
    nevNx(:, 3) = nevNx(:, 3) + nev(end, 3) + 1;
    nev = [nev; nevNx];
    waves = [waves, wavesNx];
end
% Convert Nevs to dat so that its easier to work with in Matlab
datBase = nev2dat(nevBase, 'nevreadflag', true);
% Run NASNET on NEVs
[~,nevLabelledData] = runNASNet({nev, waves},gamma, 'netFolder', netFolder, 'netname', nasNetName, 'labelSpikesAsWithWrite', true);
%% Initialize BCI-specific parameters
% Setup filepaths to where BCI decoder will be saved.
bciDecoderSaveDrive = params.bciDecoderBasePathDataComputer;
bciDecoderRelativeSaveFolder = fullfile(subject);
bciDecoderSaveFolder = fullfile(bciDecoderSaveDrive, bciDecoderRelativeSaveFolder);
success = mkdir(bciDecoderSaveFolder);
if ~success
    fprintf('\nError creating new directory for BCI parameters\n')
    fprintf('\nkeyboard here...\n')
    keyboard
end
samplingRate = params.neuralRecordingSamplingFrequencyHz; % samples/s
binSizeMs = trainParams.binSizeMs; % ms
% Parameters used to identify valid channels that will be used for BCI.
frThresh = trainParams.firingRateThreshold;
ffThresh = trainParams.fanoFactorThreshold;
coincThresh = trainParams.coincThresh;
coincTimeMs = trainParams.coincTimeMs;
timeDecoderDelay = 0; % ms;
binDecoderDelay = timeDecoderDelay/binSizeMs;
%% Preprocess and train decoder using calibration data
% remove all non-calibration trials, as they don't have correct spiking
% data and will mess up identification of internal state axi
calibrationTrialCode = codes.BACKGROUND_PROCESS_TRIAL;
calibrationTrialsInds = cellfun(@(x) any(ismember(x(:, 2), calibrationTrialCode)), {datStruct.event});

[trimmedDat, channelsKeep] = preprocessDat(datStruct, nevLabelledData, channelNumbersUse, binSizeMs, frThresh, ffThresh, coincThresh, coincTimeMs);
trimmedDat = trimmedDat(~calibrationTrialsInds);
fprintf("\n%d channels being used\n", length(channelsKeep));

% For neural signals, look at delay period right before go cue is sent
% during calibration trials
% find region that is 250ms before the binEnd code
codeForBinEnd = codes.(trainParams.binEndCode);

% Train only on trials that have matching trainingResultCodes

% Only look at the fixation periods of trainTrials; this is because the
% number of trials with fixate codes and number of fixation_off codes aren't equal. 
[fixStartIndex, ~] = cellfun(@(x) find(x(:, 2)==codeForBinStart), {trimmedDat.event}, 'uni', 0);
[fixEndIndex, ~] = cellfun(@(x) find(x(:, 2)==codeForBinEnd), {trimmedDat.event}, 'uni', 0);

% Get all spike times and spike channels from dat
spikeTimes = cellfun(@(firstSpike, spikeTimeDiffs) [firstSpike; firstSpike+cumsum(uint64(spikeTimeDiffs))], {trimmedDat.firstspike}, {trimmedDat.spiketimesdiff}, 'uni', 0);
spikeChannels = cellfun(@(spikeInfo) spikeInfo(:, 1), {trimmedDat.spikeinfo}, 'uni', 0);
%% Training the Motivation Axis Decoder
% grab the neural activity during the delay epochs of all trials that we
% wish to train on
delayValidTrials = ~cellfun('isempty', fixEndIndex) & ~cellfun('isempty', fixStartIndex) & trainTrials;
binnedSpikesDelay = cell(size(spikeTimes));
% Bin spikes involved in delay period
binnedSpikesDelay(delayValidTrials) = cellfun(...
        @(trialSpikeTimes, chanNums, fixStartEndTime)... all the inputs
        ... adding binSizeSamples/2 ensures that the last point is used
        ... while subtracting 0.1 is for a super corner case where a spike overlaps the last bin end
        ... --histcounts2 includes this in the last bin (unlike what it does for every other bin);
        ... as spike times are integers the 0.1 subtraction changes no other bins
        histcounts2(trialSpikeTimes, chanNums, (fixStartEndTime(1):binSizeSamples:fixStartEndTime(2)+binSizeSamples/2)-0.1, 1:max(channelsKeep)+1),... histcounts2 bins things in 50ms bins and groups by channel
        spikeTimes(delayValidTrials), spikeChannels(delayValidTrials), fixStartEndTimes(delayValidTrials), 'uni', 0);
% Only keep binned counts of channels that are kept
binnedSpikesDelay(delayValidTrials) = cellfun(@(bS) bS(:, channelsKeep), binnedSpikesDelay(delayValidTrials), 'uni', 0);
allBinnedCountsFixationValidTrials = binnedSpikesDelay(delayValidTrials)';

end
