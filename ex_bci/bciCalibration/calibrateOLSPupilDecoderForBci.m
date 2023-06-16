function [decoderFileLocationAndName] = calibrateOLSPupilDecoderForBci(socketsControlComm, nevFilebase, nevFilesForTrain, trainParams, subject)

global params codes

gamma = trainParams.gamma;
channelNumbersUse = trainParams.rippleChannelNumbersInBci;
netFolder = params.nasNetFolderDataComputer;
nasNetName = trainParams.nasNetwork;


% don't double read the base if it's also a train file
if any(strcmp(nevFilesForTrain, nevFilebase))
    includeBaseForTrain = true;
else
    includeBaseForTrain = false;
end
nevFilesForTrain(strcmp(nevFilesForTrain, nevFilebase)) = [];

[nevBase,waves] = readNEV(nevFilebase);
nev = nevBase;
kr = [];

% Needed to globally align times for multiple NEV files
firstTimePtShift = 0;
fnStartTimes = firstTimePtShift;

% Read in relevant Nev Files and NS5 files 
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
datBase = nev2dat(nevBase, 'nevreadflag', true);

[~,nevLabelledData] = runNASNet({nev, waves},gamma, 'netFolder', netFolder, 'netname', nasNetName, 'labelSpikesAsWithWrite', true);


bciDecoderSaveDrive = params.bciDecoderBasePathDataComputer;
bciDecoderRelativeSaveFolder = fullfile(subject);
bciDecoderSaveFolder = fullfile(bciDecoderSaveDrive, bciDecoderRelativeSaveFolder);

samplingRate = params.neuralRecordingSamplingFrequencyHz; % samples/s
binSizeMs = trainParams.binSizeMs; % ms
msPerS = 1000; % 1000 ms in a second
binSizeSamples = binSizeMs/msPerS*samplingRate;

timeDecoderDelay = 0; % ms;
binDecoderDelay = timeDecoderDelay/binSizeMs;

% Asks the control computer for the most recent decoder relative path. In
% our case, we rely on the KalmanDecoder which is trained as the first step
% of our two step decoder.
try
    sendMessageWaitAck(socketsControlComm, 'kalmanDecoderForPlane');
    intuitiveKalmanDecoderPath = receiveMessageSendAck(socketsControlComm)
catch
    warning('No decoder file was received from the control computer... Entering debugging mode')
    kalmanDecoderPath = 'satchel\Sa220922KalmanBci_12-26-46.mat';
    fprintf('Using the following decoder file: %s', kalmanDecoderPath)
    intuitiveKalmanDecoderPath = kalmanDecoderPath;
end
% Retrieve the model parameters (FA projection matrices and Kalman) from the Kalman decoder path
kalmanDecoderFilepath = fullfile(bciDecoderSaveDrive, intuitiveKalmanDecoderPath);
modelParams = load(kalmanDecoderFilepath, 'M2', 'M1', 'M0', 'K', 'beta', 'estParams', 'channelsKeep', 'zScoreLatentMat', 'zScoreSpikesMat');
channelsKeep = modelParams.channelsKeep;
channelsKeepWithDig = [0 channelsKeep];
nevLabelledData = nevLabelledData(ismember(nevLabelledData(:, 1), channelsKeepWithDig), :);

nevFilesForTrainWithBase = cat(1, {nevFilebase}, nevFilesForTrain);
trimmedDat = nev2dat(nevLabelledData, 'nevreadflag', true, 'nevfilename', nevFilesForTrainWithBase, ...
    'fnStartTimes', fnStartTimes ,'convertEyes', true, 'readNS5', true, ...
    'allowNevPause', true);
fprintf("\n%d channels being used\n", length(channelsKeep));

% NevBase is always sent; this just checks that the current file also
% includes NEVbase or not.
if ~includeBaseForTrain
    trimmedDat = trimmedDat(length(datBase)+1:end);
    nevLabelledData = nevLabelledData(size(nevBase, 1)+1:end, :);
end

% For neural signals, look at fixation period which starts when animal
% attains fixation and fix_off as when fixation dot turns off
codeForBinStart = codes.(trainParams.binStartCode);
codeForBinEnd = codes.(trainParams.binEndCode);

% Train only on trials that have matching trainingResultCodes
if iscell(trainParams.trainingResultCodes)
    resultCodesForTrialsToKeep = cellfun(@(resCode) codes.(resCode),  trainParams.trainingResultCodes);
elseif ischar(trainParams.trainingResultCodes)
    resultCodesForTrialsToKeep = codes.(trainParams.trainingResultCodes);
end
trainTrials = cellfun(@(x) any(ismember(x(:, 2), resultCodesForTrialsToKeep)), {trimmedDat.event});

% Only look at the fixation periods of trainTrials; this is because the
% number of trials with fixate codes and number of fixation_off codes aren't equal. 
[fixStartIndex, ~] = cellfun(@(x) find(x(:, 2)==codeForBinStart), {trimmedDat.event}, 'uni', 0);
[fixEndIndex, ~] = cellfun(@(x) find(x(:, 2)==codeForBinEnd), {trimmedDat.event}, 'uni', 0);

% For some reason, some trials have two fix_off codes sent, if so, select
% the first fixation code.
for i = 1:length(fixEndIndex)
    if length(fixEndIndex{i}) > 1
        fixEndIndex{i} = [fixEndIndex{i}(1)];
    end
end
% Contains lstart/end sample numbers for fixation periods
fixStartEndTimes = cellfun(@(evtArr, stInd, endInd) double(evtArr([stInd endInd], 3)), {trimmedDat.event}, fixStartIndex, fixEndIndex, 'uni', 0);

% Get all spike times and spike channels from dat
spikeTimes = cellfun(@(firstSpike, spikeTimeDiffs) [firstSpike; firstSpike+cumsum(uint64(spikeTimeDiffs))], {trimmedDat.firstspike}, {trimmedDat.spiketimesdiff}, 'uni', 0);
spikeChannels = cellfun(@(spikeInfo) spikeInfo(:, 1), {trimmedDat.spikeinfo}, 'uni', 0);
%% Training the Pupil Axis Decoder
% grab the neural activity during the fixation periods of all calibration
% trials that we want to train on
fixValidTrials = ~cellfun('isempty', fixEndIndex) & ~cellfun('isempty', fixStartIndex) & trainTrials;
binnedSpikesFixation = cell(size(spikeTimes));
% Bin spikes involved in fixation period
binnedSpikesFixation(fixValidTrials) = cellfun(...
        @(trialSpikeTimes, chanNums, fixStartEndTime)... all the inputs
        ... adding binSizeSamples/2 ensures that the last point is used
        ... while subtracting 0.1 is for a super corner case where a spike overlaps the last bin end
        ... --histcounts2 includes this in the last bin (unlike what it does for every other bin);
        ... as spike times are integers the 0.1 subtraction changes no other bins
        histcounts2(trialSpikeTimes, chanNums, (fixStartEndTime(1):binSizeSamples:fixStartEndTime(2)+binSizeSamples/2)-0.1, 1:max(channelsKeep)+1),... histcounts2 bins things in 50ms bins and groups by channel
        spikeTimes(fixValidTrials), spikeChannels(fixValidTrials), fixStartEndTimes(fixValidTrials), 'uni', 0);
% Only keep binned counts of channels that are kept
binnedSpikesFixation(fixValidTrials) = cellfun(@(bS) bS(:, channelsKeep), binnedSpikesFixation(fixValidTrials), 'uni', 0);
allBinnedCountsFixationValidTrials = binnedSpikesFixation(fixValidTrials)';

% Grab the pupil activity during fixation
binnedPupilSizesFixation = cell(size(spikeTimes));
binnedPupilSizesFixation(fixValidTrials) = cellfun(...
    @(trialPupilData,fixStartEndTime)
    fixStartEndTime(
%%
% output: directions for task in save file
subjectCamelCase = lower(subject);
subjectCamelCase(1) = upper(subjectCamelCase(1));
bciDecoderSaveName = sprintf('%s%sNeuralEngagementBci_%s.mat', subjectCamelCase(1:2), datestr(today, 'yymmdd'), datestr(now, 'HH-MM-SS'));
save(fullfile(bciDecoderSaveFolder, bciDecoderSaveName), 'M0', 'M1', 'M2','interpolatedConditionMeans', 'interpolatedNeuralEngagementAxes' ,'interpolatedNeuralEngagementValueMeans','interpolatedNeuralEngagementValueStds', 'zScoreSpikesMat','d','betaOrthNorm', 'channelsKeep', 'nevFilebase', 'nevFilesForTrain', 'includeBaseForTrain', 'nasNetName', 'intuitiveKalmanDecoderPath');
decoderFileLocationAndName = fullfile(bciDecoderRelativeSaveFolder, bciDecoderSaveName);
%% Helper function for extracting pupil data for each trial
function binTrialPupilData(trialPupilData, fixationStartEndTimes)
    trialStartSample = trialPupilData.startsample;
    trialPupilSizes = trialPupilData.trial;
    % Convert to indices by subtracting startsample from startEndTimes
    trialFixationStartEndIndices = fixationStartEndTimes - trialStartSample;
    % Extract pupil size samples that fall within startEndTimes
    pupilSizesDuringPeriod = trialPupilSizes(trialFixationStartEndIndices(1):trialFixationStartEndIndices(2));
    binStartEndIndices = (trialFixationStartEndIndices(1): binSizeSamples: trialFixationStartEndIndices(2) + binSizeSamples/2) - 0.1;
end

