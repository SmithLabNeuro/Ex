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
% Base usually set to the first file of training Nevs
% TODO: Figure out what this nevFileBase actually does ( Maybe used to
% determine timepoints w.r.t other Nev files)
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
% Run NASNET on NEVs
[~,nevLabelledData] = runNASNet({nev, waves},gamma, 'netFolder', netFolder, 'netname', nasNetName, 'labelSpikesAsWithWrite', true);
% Convert Nevs to dat so that its easier to work with in Matlab
datStruct = nev2dat(nevLabelledData, 'nevreadflag', true);
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
% TODO: Figure these out
% Parameters used to identify valid channels that will be used for BCI.
frThresh = trainParams.firingRateThreshold;
ffThresh = trainParams.fanoFactorThreshold;
coincThresh = trainParams.coincThresh;
coincTimeMs = trainParams.coincTimeMs;
%% Preprocess and train decoder using calibration data
% remove all non-calibration trials, as they don't have correct spiking
% data and will mess up identification of internal state axis
nonCalibrationTrialCode = codes.BACKGROUND_PROCESS_TRIAL; % Code 250 should refer to nonCalibration Trials
nonCalibrationTrialsInds = cellfun(@(x) any(ismember(x(:, 2), nonCalibrationTrialCode)), {datStruct.event});

% Preprocess data according to criteria we apply (FR threshold, coincidence
% detection, etc.)
[trimmedDat, channelsKeep] = preprocessDat(datStruct, nevLabelledData, channelNumbersUse, binSizeMs, frThresh, ffThresh, coincThresh, coincTimeMs);
% Dat should only contain calibration trials
trimmedDat = trimmedDat(~nonCalibrationTrialsInds);
fprintf("\n%d channels being used\n", length(channelsKeep));

% For neural signals, look at delay period right before go cue is sent
% during calibration trials
% TODO: Specify what this calBinEndCode will be for our task; it's currently
% used here to determine the end of our bins that we'll use for training
codeForBinEnd = codes.(trainParams.calBinEndCode); 

% Train only on trials that have matching trainingResultCodes
% TODO: Specify what trainingResultCodes we will use in our task
if iscell(trainParams.trainingResultCodes)
    resultCodesForTrialsToKeep = cellfun(@(resCode) codes.(resCode),  trainParams.trainingResultCodes);
elseif ischar(trainParams.trainingResultCodes)
    resultCodesForTrialsToKeep = codes.(trainParams.trainingResultCodes);
end

% Identify calibration trial indices that have matching result codes
trainTrials = cellfun(@(x) any(ismember(x(:, 2), resultCodesForTrialsToKeep)), {trimmedDat.event});

% Only look at the delay periods of trainTrials; this will be timeBeforeEndCode ms before
% the calBinEndCode that we provide
% Identify samples that correspond to the end and beginning for epoch
% Go back 250ms from end Sample 
% TODO: Is this the time period we want before the epoch?
timeFromEndOfEpoch = 250; %in ms
samplingRate = params.neuralRecordingSamplingFrequencyHz; % samples/s
% First element of vector in cellfun is the start sample for delay epoch
% and second sample is the end of that epoch
[epochStartEndSampleIndices] = cellfun(@(x) [x(x(:, 2)==codeForBinEnd, 3) - samplingRate*timeFromEndOfEpoch/1000, x(x(:, 2)==codeForBinEnd, 3)] , {trimmedDat.event}, 'uni', 0);
% Get all spike times and spike channels from dat
spikeTimes = cellfun(@(firstSpike, spikeTimeDiffs) [firstSpike; firstSpike+cumsum(uint64(spikeTimeDiffs))], {trimmedDat.firstspike}, {trimmedDat.spiketimesdiff}, 'uni', 0);
spikeChannels = cellfun(@(spikeInfo) spikeInfo(:, 1), {trimmedDat.spikeinfo}, 'uni', 0);
%% Grab Neural activity that will be used to fit the internal state axis 
% grab the neural activity during the delay epochs for training trials
% All trials that have an epoch end are valid trials to train on
binSizeMs = trainParams.binSizeMs; % ms
msPerS = 1000; % 1000 ms in a second
binSizeSamples = binSizeMs/msPerS*samplingRate;
delayValidTrials = ~cellfun('isempty', epochStartEndSampleIndices) & trainTrials;
% Initialize cell that will store delay epoch spikes
binnedSpikesDelay = cell(size(spikeTimes));
% Bin spikes involved in delay period
binnedSpikesDelay(delayValidTrials) = cellfun(...
        @(trialSpikeTimes, chanNums, trialEpochStartEndTime)... all the inputs
        ... adding binSizeSamples/2 ensures that the last point is used
        ... while subtracting 0.1 is for a super corner case where a spike overlaps the last bin end
        ... --histcounts2 includes this in the last bin (unlike what it does for every other bin);
        ... as spike times are integers the 0.1 subtraction changes no other bins
        histcounts2(trialSpikeTimes, chanNums, (trialEpochStartEndTime(1):binSizeSamples:trialEpochStartEndTime(2)+binSizeSamples/2)-0.1, 1:max(channelsKeep)+1),... histcounts2 bins things in 50ms bins and groups by channel
        spikeTimes(delayValidTrials), spikeChannels(delayValidTrials), epochStartEndSampleIndices(delayValidTrials), 'uni', 0);
% Only keep binned counts of channels that are kept
binnedSpikesDelay(delayValidTrials) = cellfun(@(bS) bS(:, channelsKeep), binnedSpikesDelay(delayValidTrials), 'uni', 0);
allTrialsDelayEpochBinnedCounts = binnedSpikesDelay(delayValidTrials)';
% TODO: Do we take the average of the bins for the actual BCI? How do we use these bins?
%% Fit FA to binned counts to help denoise 
% Concatenate all trials binned spikes
binnedSpikesAllConcat = cat(1, allTrialsDelayEpochBinnedCounts{:});
% TODO: Figure out the number of latents necessary
numLatents = trainParams.numberFaLatents; % this is how many latents we'll project to, no matter what...
if trainParams.zScoreSpikes
    zScoreSpikesMat = diag(1./std(binnedSpikesAllConcat, [], 1));
    % Subtract the mean spike count vector and then divide by the variance
    % of spikes
    binnedSpikesAllConcat = (binnedSpikesAllConcat - mean(binnedSpikesAllConcat)) * zScoreSpikesMat;
else
    % Do nothing
end
% Identify FA Space using ALL bins of delay period 
% (num_bins_in_delay* num_trials x num_channels)
[estParams, ~] = fastfa(binnedSpikesAllConcat', numLatents);
% Beta is simply the projection matrix into the factor space
[~, ~, beta] = fastfa_estep(binnedSpikesAllConcat', estParams);
%% Fit LDA to these FA Projections
% Keep trial parameters of valid trials that will be used for training LDA
delayValidTrialParams = trimmedDat(delayValidTrials);
numValidTrials = length(delayValidTrialParams);
validTrialRewardLabels = nan(numValidTrials, 1);
for k = 1:numValidTrials
    validTrialRewardLabels(k) = delayValidTrialParams(k).params.trial.variableRewardIdx;
end
% Take Mean binned activity during delay period then project into FA space
meanDelayValidTrialFAProjs = cellfun(@(x) beta*(mean(x)' - estParams.d), allTrialsDelayEpochBinnedCounts, 'UniformOutput', false);
% Concatenate across all trials to get a num_valid_trials x num_fa_latents
ldaTrainX = cat(2, meanDelayValidTrialFAProjs{:})';
% Fit LDA on FA projections
ldaParams = fit_LDA(ldaTrainX, validTrialRewardLabels);
%% Save model parameters 
bciDecoderSaveName = sprintf('rewardAxisDecoder_%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
save(fullfile(bciDecoderSaveFolder, bciDecoderSaveName), 'ldaParams', 'beta','channelsKeep', 'nevFilebase', 'nevFilesForTrain', 'includeBaseForTrain', 'nasNetName');
decoderFileLocationAndName = fullfile(bciDecoderRelativeSaveFolder, bciDecoderSaveName);
fprintf('decoder file saved at : %s\n', decoderFileLocationAndName)
end
