function [decoderFileLocationAndName] = calibrateMultipleAxisDecoderForBci(~, nevFilebase, nevFilesForTrain, trainParams, subject, offlineFlag)

% inputs:
% trainParams - comes from bci_rewardAxisDecoder.xml

global params codes
% NASNet variables
gamma = trainParams.gamma;
nasNetName = trainParams.nasNetwork;
% If you ever see an empty params, make sure to run recordEX first
netFolder = params.nasNetFolderDataComputer;
%netFolder = 'C:\Users\rigmdata\Documents\spikesort\nasnet\networks';
% Specify channels that will be used for BCI
channelNumbersUse = trainParams.rippleChannelNumbersInBci;

% don't double read the base if it's also a train file
if any(strcmp(nevFilesForTrain, nevFilebase))
    includeBaseForTrain = true;
else
    includeBaseForTrain = false;
end
nevFilesForTrain(strcmp(nevFilesForTrain, nevFilebase)) = [];

if nargin < 6
    offlineFlag = false;
end

%% Read in NEV files that will be used for training our decoder
% Base usually set to the first file of training Nevs
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
codeForBinStart = codes.(trainParams.calBinStartCode);
% Keep track of how many bins after start of epoch we want to actually use
timeFromStartOfEpoch = trainParams.timeFromStartOfEpoch;
codeForBinEnd = codes.(trainParams.calBinEndCode); 

% Train only on trials that have matching trainingResultCodes
if iscell(trainParams.trainingResultCodes)
    resultCodesForTrialsToKeep = cellfun(@(resCode) codes.(resCode),  trainParams.trainingResultCodes);
elseif ischar(trainParams.trainingResultCodes)
    resultCodesForTrialsToKeep = codes.(trainParams.trainingResultCodes);
end

% Identify calibration trial indices that have matching result codes
trainTrials = cellfun(@(x) any(ismember(x(:, 2), resultCodesForTrialsToKeep)), {trimmedDat.event});

% Identify samples that correspond to the end and beginning for epoch
samplingRate = params.neuralRecordingSamplingFrequencyHz; % samples/s
% First element of vector in cellfun is the start sample for delay epoch
% and second sample is the end of that epoch
[epochStartEndSampleIndices] = cellfun(@(x) [x(x(:, 2)==codeForBinStart, 3) + samplingRate*timeFromStartOfEpoch/1000, x(x(:, 2)==codeForBinEnd, 3)] , {trimmedDat.event}, 'uni', 0);
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

delayValidTrialParams = trimmedDat(delayValidTrials);
numValidTrials = length(delayValidTrialParams);
validTrialRewardLabels = nan(numValidTrials, 1);
validTrialTargetLabels = nan(numValidTrials, 1);
for k = 1:numValidTrials
    validTrialRewardLabels(k) = delayValidTrialParams(k).params.trial.variableRewardIdx;
    validTrialTargetLabels(k) = delayValidTrialParams(k).params.trial.targetAngle;
end
%% Fit FA to binned counts to help denoise 
% Concatenate all trials binned spikes
binnedSpikesAllConcat = cat(1, allTrialsDelayEpochBinnedCounts{:});
numLatents = trainParams.numberFaLatents; % this is how many latents we'll project to, no matter what...
% For this task, we should NOT use Z-scoring
if trainParams.zScoreSpikes
    % Corresponds to 1/sigma for each neurons
    zScoreSpikesMat = diag(1./std(binnedSpikesAllConcat, [], 1));
    % Corresponds to mu/sigma for each neuron (should be a 1 x
    % num_channels)
    zScoreSpikesMuTerm =  mean(binnedSpikesAllConcat)./std(binnedSpikesAllConcat, [], 1);
else
    % Identity matrix if not z-scoring
    zScoreSpikesMat =  eye(size(binnedSpikesAllConcat, 2));
    % zeros if not z-scoring
    zScoreSpikesMuTerm = zeros(size(mean(binnedSpikesAllConcat))); 
end
% Subtract the mean spike count vector and then divide by the variance
% of spikes
binnedSpikesAllConcat = binnedSpikesAllConcat*zScoreSpikesMat - zScoreSpikesMuTerm;
zScoreSpikesMuTerm = zScoreSpikesMuTerm'; % transpose for decoder
% Identify FA Space using ALL bins of delay period
% (num_bins_in_delay* num_trials x num_channels)
[estFAParams, ~] = fastfa(binnedSpikesAllConcat', numLatents);
% Beta is simply the projection matrix into the factor space
[~, ~, beta] = fastfa_estep(binnedSpikesAllConcat', estFAParams);
%% Temporally Smooth Fa Projections
alpha = trainParams.alpha;
% Find FA projections for all trials
faProjsByTrial = cellfun(@(x) beta*(x' - estFAParams.d), allTrialsDelayEpochBinnedCounts, 'UniformOutput', false);
% Exponentially smooth each trial's bins, include previous bins' effects.
% Seed this way for smoothing the FA latents
initialSeedValue = mean(horzcat(faProjsByTrial{:}),2); % Should be zero vector
smoothedFaProjsByTrial = cellfun(@(x) exponentialSmoother(x, alpha, initialSeedValue), faProjsByTrial, 'UniformOutput', false);
validTrialRewardRepeatedLabels = [];
validTrialTargRepeatedLabels = [];
for k = 1:length(smoothedFaProjsByTrial)
    numBinsInTrial = size(smoothedFaProjsByTrial{k},2);
    validTrialTargRepeatedLabels = [validTrialTargRepeatedLabels,  repmat(validTrialTargetLabels(k), 1,numBinsInTrial)];
    validTrialRewardRepeatedLabels = [validTrialRewardRepeatedLabels,  repmat(validTrialRewardLabels(k), 1,numBinsInTrial)];
end
%% Fit LDA to smoothed FA Projections
% Concatenate across all trials to get a num_valid_trials x num_fa_latents
ldaTrainX = cat(2, smoothedFaProjsByTrial{:})';
% Fit LDA on FA projections to find reward axis
rewardLdaParams = fit_LDA(ldaTrainX, validTrialRewardRepeatedLabels);
%% Get reward axis means and magnitude
% Reorient projection vector such that the largest reward value is the
% largest value on this axis
% Keep track of Small/large 
smallRewardIndices = find(validTrialRewardRepeatedLabels == 1);
largeRewardIndices = find(validTrialRewardRepeatedLabels == 3);
disp ('Checking Reward axis for flip')

[~, ~, ~, ~, ~,~, rewardLDAParams] = flipAxesBasedOnCondition(smallRewardIndices, largeRewardIndices, rewardLdaParams);
%% Fit LDA to find axes based on smoothed FA Projections
% target axis
% 1) subsample trials to get only up and down targets
%targetsToSeparate = [90, 270; 0, 180]; % num_axes x num_targets
targetsToSeparate = [90,270];
multipleAxesLDAParams = [];
mappingTargetParams = {};
% Store target specific parameters that will be used in the decoder
for k = 1:size(targetsToSeparate,1)
    trialsWithFirstTargIdx = find(validTrialTargetLabels == targetsToSeparate(k,1));
    trialsWithSecondTargIdx = find(validTrialTargetLabels == targetsToSeparate(k,2));
    faProjsWithFirstTarg = cell2mat({smoothedFaProjsByTrial{trialsWithFirstTargIdx}});
    faProjsWithSecondTarg = cell2mat({smoothedFaProjsByTrial{trialsWithSecondTargIdx}});
    % Setting label for first target angle to be 1 and second to 2
    firstAngleLabels = repmat(targetsToSeparate(k,1), 1,size(faProjsWithFirstTarg,2));
    secondAngleLabels = repmat(targetsToSeparate(k,2), 1,size(faProjsWithSecondTarg,2));
    % Stack our training samples
    trainX = [faProjsWithFirstTarg, faProjsWithSecondTarg]';
    trainY = [firstAngleLabels, secondAngleLabels];
    % Keep track of LDA params
    axisLdaParams = fit_LDA(trainX, trainY);
    firstTargIndices = find(trainY == targetsToSeparate(k,1));
    secondTargIndices = find(trainY == targetsToSeparate(k,2));
    fprintf('Current axis for separating %i and %i\n', targetsToSeparate(k,1), targetsToSeparate(k,2))
    % Flip if necessary
    [firstTargProjs, secondTargProjs, firstTargMean, secondTargMean, firstTargProjsSD,secondTargProjsSD, axisLdaParams] = flipAxesBasedOnCondition(firstTargIndices, secondTargIndices, axisLdaParams);
    targetSpecificParams = struct('angle', {}, 'mean', {}, 'std', {}, 'LDAProjs', {});
    % Add struct array for first target
    targetSpecificParams(end + 1) = struct('angle', targetsToSeparate(k,1), 'mean', firstTargMean, 'std',firstTargProjsSD, 'LDAProjs',  firstTargProjs);
    % Add struct array for second target
    targetSpecificParams(end + 1) = struct('angle', targetsToSeparate(k,2), 'mean', secondTargMean, 'std',secondTargProjsSD, 'LDAProjs',  secondTargProjs);
    % Plot projections along this 1D intuitive axis
    figure;
    hold on
    histogram(firstTargProjs, 'binWidth', 0.25, 'FaceColor', 'b' ,'DisplayName', 'Targ1')
    histogram(secondTargProjs, 'binWidth', 0.25, 'FaceColor', 'r' ,'DisplayName', 'Targ2')
    xline(firstTargMean, '--b', 'LineWidth', 2, 'DisplayName', 'Mean of Targ 1 Projs')
    xline(secondTargMean, '--r', 'LineWidth', 2, 'DisplayName', 'Mean of Targ 2 Projs')
    hold off
    title(sprintf('Projections of calibration bins along 1D intuitive where Targ1: %i Targ2: %i', targetsToSeparate(k,1), targetsToSeparate(k,2)))
    legend()
    % Identify ranges that will be used for each state 
    multipleAxesLDAParams = [multipleAxesLDAParams, axisLdaParams];
    mappingTargetParams{end+1} = targetSpecificParams;
end
% Decoder should track: targetStates, multipleAxesLDAParams, and the ranges
%% Save model parameters 
subjectCamelCase = lower(subject);
subjectCamelCase(1) = upper(subjectCamelCase(1));
% bciDecoderSaveName = sprintf('rewardAxisDecoder_%s.mat', datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
if offlineFlag
    % Use provided nev filename for date
    fileBackSlashIndices = find(nevFilebase == '\');
    if length(fileBackSlashIndices) > 1
        fileBackSlashIndex = fileBackSlashIndices(end);
    else
        fileBackSlashIndex = fileBackSlashIndices;
    end
    bciDecoderSaveName = sprintf('%sIntuitiveAxesDecoderBci_%s_offline.mat',nevFilebase(fileBackSlashIndex+1:end-4) , datestr(now, 'dd-mmm-yyyy_HH-MM-SS') );
else
    % Use current date as naming convention for online
    bciDecoderSaveName = sprintf('%s%sIntuitiveAxesDecoderBci_%s.mat', subjectCamelCase(1:2), datestr(today, 'yymmdd'), datestr(now, 'HH-MM-SS'));
end
% Keep track of reward Axes activity just to see how neural activity acts
% offline
save(fullfile(bciDecoderSaveFolder, bciDecoderSaveName), ...
'estFAParams', 'beta', 'zScoreSpikesMat', 'zScoreSpikesMuTerm', 'channelsKeep', ...
'nevFilebase', 'nevFilesForTrain', 'includeBaseForTrain', 'nasNetName', 'rewardLDAParams', ...
'mappingTargetParams', 'multipleAxesLDAParams' , 'trainParams', 'delayValidTrials', 'binnedSpikesDelay');

decoderFileLocationAndName = fullfile(bciDecoderRelativeSaveFolder, bciDecoderSaveName);
fprintf('decoder file saved at : %s\n', decoderFileLocationAndName)
end