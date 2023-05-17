function [outputStruct] = trainRewardAxisDecoderFromNev(nevFilebase, subject)
% function [M0, M1, M2, channelsKeep, A, Q, C, R, beta, K] = trainKalmanDecoderFromNev(nev, trainParams, params, codes, datBase, nevBase,...
%         netLabels, gamma, channelNumbersUse, Qvalue, numLatents, velocityToCalibrateWith, includeBaseForTrain)

exGlobals;
dataPath = 'E:\';
xmlFile = 'bci_rewardAxisDecoder.xml';
[~,machineInit] = system('hostname');
machine = lower(deblank(cell2mat(regexp(machineInit, '^[^\.]+', 'match'))));

nevFilename = fullfile(dataPath, nevFilebase);
[~, trainParams, ~, ~] = readExperiment(xmlFile,subject,machine);

% NASNet variables
gamma = trainParams.gamma;
nasNetName = trainParams.nasNetwork;
netFolder = params.nasNetFolderDataComputer;
% Specify channels that will be used for BCI
channelNumbersUse = trainParams.rippleChannelNumbersInBci;

%% Read in NEV files that will be used for training our decoder
[nev,waves] = readNEV(nevFilename);

% Run NASNET on NEVs
[~,nevLabelledData] = runNASNet({nev, waves}, gamma, 'netFolder', netFolder, 'netname', nasNetName, 'labelSpikesAsWithWrite', true);
% Convert Nevs to dat so that its easier to work with in Matlab
datStruct = nev2dat(nevLabelledData, 'nevreadflag', true);
%% Initialize BCI-specific parameters
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
timeFromEndOfEpoch = trainParams.timeBeforeEndCode; %in ms
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

% subsample trials to get same number of rewards
delayValidTrialParams = trimmedDat(delayValidTrials);
numValidTrials = length(delayValidTrialParams);
validTrialRewardLabels = nan(numValidTrials, 1);
for k = 1:numValidTrials
    validTrialRewardLabels(k) = delayValidTrialParams(k).params.trial.variableRewardIdx;
end
uniqueRewardIdxs = unique(validTrialRewardLabels);
numTrialsPerReward = histc(validTrialRewardLabels, uniqueRewardIdxs);
minNumTrialsPerReward = min(numTrialsPerReward);
subsampleIdx = zeros(numValidTrials,1);
for i=1:length(uniqueRewardIdxs)
    r = uniqueRewardIdxs(i);
    n = numTrialsPerReward(i);
    fullRewardIdx = validTrialRewardLabels == r;
    if n > minNumTrialsPerReward
        trimIdx = find(cumsum(fullRewardIdx)==minNumTrialsPerReward+1);
        fullRewardIdx(trimIdx:end) = false;
    end
    subsampleIdx = subsampleIdx | fullRewardIdx;
end
% sub sample
allTrialsDelayEpochBinnedCounts = allTrialsDelayEpochBinnedCounts(subsampleIdx);
validTrialRewardLabels = validTrialRewardLabels(subsampleIdx);
%% Fit FA to binned counts to help denoise 
% Concatenate all trials binned spikes
binnedSpikesAllConcat = cat(1, allTrialsDelayEpochBinnedCounts{:});
% TODO: Figure out the number of latents necessary
numLatents = trainParams.numberFaLatents; % this is how many latents we'll project to, no matter what...
if trainParams.zScoreSpikes
    % Corresponds to 1/sigma for each neurons
    zScoreSpikesMat = diag(1./std(binnedSpikesAllConcat, [], 1));
    % Corresponds to mu/sigma for each neuron
    zScoreSpikesMuTerm =  mean(binnedSpikesAllConcat)./std(binnedSpikesAllConcat, [], 1);
else
    % Identity matrix if not z-scoring
    zScoreSpikesMat =  eye(size(binnedSpikesAllConcat, 2));
    % zeros if not z-scoring
    zScoreSpikesMuTerm = zeros(mean(binnedSpikesAllConcat)); 
end
% Subtract the mean spike count vector and then divide by the variance
% of spikes
binnedSpikesAllConcat = binnedSpikesAllConcat*zScoreSpikesMat - zScoreSpikesMuTerm;
% Identify FA Space using ALL bins of delay period
% (num_bins_in_delay* num_trials x num_channels)
[estFAParams, ll] = fastfa(binnedSpikesAllConcat', numLatents);
% Beta is simply the projection matrix into the factor space
[~, ~, beta] = fastfa_estep(binnedSpikesAllConcat', estFAParams);
%% Temporally Smooth Fa Projections
alpha = trainParams.alpha;
numBinsForLDA = trainParams.numCalBins; % Specify number of calibration bins that will be used for training LDA
% Find FA projections for all trials
faProjsByTrial = cellfun(@(x) beta*(x' - estFAParams.d), allTrialsDelayEpochBinnedCounts, 'UniformOutput', false);
% Exponentially smooth each trial's bins, include previous bins' effects 
smoothedFaProjsByTrial = cellfun(@(x) exponentialSmoother(x, alpha), faProjsByTrial, 'UniformOutput', false);
% Select the last numBinsForLDA for each calibration trial
smoothedFaProjsByTrial = cellfun(@(x) x(:, end-(numBinsForLDA-1):end), smoothedFaProjsByTrial, 'UniformOutput', false);
% Repeat reward labels so that each bin in a trial has a label for LDA
% Repmat replicates by column and reshape reorders row-wise
binRewardLabelsByTrial = reshape(repmat(validTrialRewardLabels,numBinsForLDA,1),length(validTrialRewardLabels)*numBinsForLDA ,1);
%% Fit LDA to smoothed FA Projections
% Concatenate across all trials to get a num_valid_trials x num_fa_latents
ldaTrainX = cat(2, smoothedFaProjsByTrial{:})';
% Fit LDA on FA projections
ldaParams = fit_LDA(ldaTrainX, binRewardLabelsByTrial);
%% Get reward axis means and magnitude
% Reorient projection vector such that the largest reward value is the
% largest value on this axis
largeRewardTrialIndices = validTrialRewardLabels(validTrialRewardLabels == 3);
smallRewardTrialIndices = validTrialRewardLabels(validTrialRewardLabels == 1);
largeRewardMeanProj = mean(ldaParams.projData(largeRewardTrialIndices));
smallRewardMeanProj = mean(ldaParams.projData(smallRewardTrialIndices));
% Flip reward axis only if smallReward projection is higher than large reward projection
if smallRewardMeanProj > largeRewardMeanProj
    ldaParams.projVec = ldaParams.projVec*-1;
    ldaParams.projData = ldaParams.projData*-1;
    largeRewardMeanProj = -largeRewardMeanProj;
    smallRewardMeanProj = -smallRewardMeanProj;
end

rewardAxisRange = largeRewardMeanProj - smallRewardMeanProj; % used to calculate current neural distance

%% Save model parameters 
outputStruct = struct();
outputStruct.ldaParams = ldaParams;
outputStruct.estFAParams = estFAParams;
outputStruct.beta = beta;
outputStruct.zScoreSpikesMat = zScoreSpikesMat;
outputStruct.zScoreSpikesMuTerm = zScoreSpikesMuTerm;
outputStruct.channelsKeep = channelsKeep;
outputStruct.largeRewardMeanProj = largeRewardMeanProj;
outputStruct.smallRewardMeanProj = smallRewardMeanProj;
outputStruct.rewardAxisRange = rewardAxisRange;

end
function smoothedBins = exponentialSmoother(currTrialBins, alpha)
% Assumes currTrialBins is of num_dims x num_time_bins
    smoothedBins = nan(size(currTrialBins));
    for i =1:size(smoothedBins,2)
        % Set first value to currTrialBins first value
        if i == 1
            smoothedBins(:,i) = currTrialBins(:,i);
        else
            smoothedBins(:,i) = (1-alpha)*smoothedBins(:,i-1) + alpha*currTrialBins(:,i);
        end
    end
end