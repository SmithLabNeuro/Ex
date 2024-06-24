function decoderFileLocationAndName = calibrateAxisForClosedLoop(socketsControlComm, nevFilebase, nevFilesForTrain, trainParams, subject)

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


% Load nev file and create dat structure
[nevBase,waves] = readNEV(nevFilebase);
nev = nevBase;
for nevFlInd = 1:length(nevFilesForTrain)
    [nevNx, wavesNx] = readNEV(nevFilesForTrain{nevFlInd});
    nevNx(:, 3) = nevNx(:, 3) + nev(end, 3) + 1;
    nev = [nev; nevNx];
    waves = [waves, wavesNx];
end
datBase = nev2dat(nevBase, 'nevreadflag', true);

[slabel,nevLabelledData] = runNASNet({nev, waves},gamma, 'netFolder', netFolder, 'netname', nasNetName);
datStruct = nev2dat(nevLabelledData, 'nevreadflag', true);

if ~includeBaseForTrain
    datStruct = datStruct(length(datBase)+1:end);
    nevLabelledData = nevLabelledData(size(nevBase, 1)+1:end, :);
end


% Prepare a folder to store the calibration parameter file
bciDecoderSaveDrive = params.bciDecoderBasePathDataComputer;
bciDecoderRelativeSaveFolder = fullfile(subject);
bciDecoderSaveFolder = fullfile(bciDecoderSaveDrive, bciDecoderRelativeSaveFolder);
success = mkdir(bciDecoderSaveFolder);
if ~success
    fprintf('\nError creating new directory for BCI parameters\n')
    fprintf('\nkeyboard here...\n')
    keyboard
end


% Define calibration related parameters
frThresh = trainParams.firingRateThreshold;
ffThresh = trainParams.fanoFactorThreshold;
coincThresh = trainParams.coincThresh;
coincTimeMs = trainParams.coincTimeMs;

samplingRate = params.neuralRecordingSamplingFrequencyHz; % samples/s
binSizeMs = trainParams.binSizeMs; % ms
msPerS = 1000; % 1000 ms in a second
binSizeSamples = binSizeMs/msPerS*samplingRate;

% this is the code after which the neural activity data is sent to bci computer; 
% could be FIX_OFF, or TARG_OFF, or any of the existing codes
neuralSignalForBciStart = codes.(trainParams.calibrateStartCode);% [CURSOR_ON] expParams.bciStartCode;
% neural signal no longer controls BCI
neuralSignalForBciEnd = codes.(trainParams.calibrateEndCode);% [CURSOR_OFF CORRECT] expParams.bciEndCodes;

% these are the codes to indicate the stard/end of pre uStim period
neuralSignalForPreStart = codes.(trainParams.preActStartCode);
neuralSignalForPreEnd = codes.(trainParams.preActEndCode);

% these are the codes to indicate the stard/end of pre uStim period
neuralSignalForPostStart = codes.(trainParams.postActStartCode);
neuralSignalForPostEnd = codes.(trainParams.postActEndCode);

% Run preprocess
[trimmedDat, channelsKeep] = preprocessDatWithNoStimTrial(datStruct, nevLabelledData, channelNumbersUse, binSizeMs, frThresh, ffThresh, coincThresh, coincTimeMs);
% [trimmedDat, channelsKeep] = preprocessDat(datStruct, nevLabelledData, channelNumbersUse, binSizeMs, frThresh, ffThresh, coincThresh, coincTimeMs);
spikeTimes = cellfun(@(firstSpike, spikeTimeDiffs) [firstSpike; firstSpike+cumsum(uint64(spikeTimeDiffs))], {trimmedDat.firstspike}, {trimmedDat.spiketimesdiff}, 'uni', 0);
spikeChannels = cellfun(@(spikeInfo) spikeInfo(:, 1), {trimmedDat.spikeinfo}, 'uni', 0);


% Define the indices and times to filter data during a specific peirod (such as pre uStim period)
[bciStartIndex, ~] = cellfun(@(x) find(x(:, 2)==neuralSignalForBciStart), {trimmedDat.event}, 'uni', 0);
[bciEndIndexInit, ~] = cellfun(@(x) find(x(:, 2)==neuralSignalForBciEnd), {trimmedDat.event}, 'uni', 0);
bciEndIndex = cellfun(@(x) min(x), bciEndIndexInit, 'uni', 0);
bciStartEndTimes = cellfun(@(evtArr, stInd, endInd) double(evtArr([stInd endInd], 3)), {trimmedDat.event}, bciStartIndex, bciEndIndex, 'uni', 0);

[preStartIndex, ~] = cellfun(@(x) find(x(:, 2)==neuralSignalForPreStart), {trimmedDat.event}, 'uni', 0);
[preEndIndexInit, ~] = cellfun(@(x) find(x(:, 2)==neuralSignalForPreEnd), {trimmedDat.event}, 'uni', 0);
preEndIndex = cellfun(@(x) min(x), preEndIndexInit, 'uni', 0);
preStartEndTimes = cellfun(@(evtArr, stInd, endInd) double(evtArr([stInd endInd], 3)), {trimmedDat.event}, preStartIndex, preEndIndex, 'uni', 0);

[postStartIndex, ~] = cellfun(@(x) find(x(:, 2)==neuralSignalForPostStart), {trimmedDat.event}, 'uni', 0);
[postEndIndexInit, ~] = cellfun(@(x) find(x(:, 2)==neuralSignalForPostEnd), {trimmedDat.event}, 'uni', 0);
postEndIndex = cellfun(@(x) min(x), postEndIndexInit, 'uni', 0);
postStartEndTimes = cellfun(@(evtArr, stInd, endInd) double(evtArr([stInd endInd], 3)), {trimmedDat.event}, postStartIndex, postEndIndex, 'uni', 0);


% Extract trial indices data (bool) used in calibration
if iscell(trainParams.trainingResultCodes)
    resultCodesForTrialsToKeep = cellfun(@(resCode) codes.(resCode),  trainParams.trainingResultCodes);
elseif ischar(trainParams.trainingResultCodes)
    resultCodesForTrialsToKeep = codes.(trainParams.trainingResultCodes);
end
trainTrialsResultCode = cellfun(@(x) any(ismember(x(:, 2), resultCodesForTrialsToKeep)), {trimmedDat.event});
trainTrialsStimChan = cellfun(@(x) any(ismember(x.trial.xippmexStimChan, trainParams.traininguStimChan)), {trimmedDat.params});
trainTrials = trainTrialsResultCode & trainTrialsStimChan;


%% Fit FA
% find FA axis along which the closed loop algorothm controls the activity

% grab the neural activity during the BCI period (i.e. from STIM2_ON to STIM2_OFF)
bciValidTrials = ~cellfun('isempty', bciEndIndex) & ~cellfun('isempty', bciStartIndex) & trainTrials;
binnedSpikesBci = cell(size(spikeTimes));
binnedSpikesBci(bciValidTrials) = cellfun(...
        @(trialSpikeTimes, chanNums, bciStartEndTime)... all the inputs
        ... adding binSizeSamples/2 ensures that the last point is used
        histcounts2(trialSpikeTimes, chanNums, bciStartEndTime(1):binSizeSamples:bciStartEndTime(2)+binSizeSamples/2, 1:max(channelsKeep)+1),... histcounts2 bins things in 50ms bins and groups by channel
        spikeTimes(bciValidTrials), spikeChannels(bciValidTrials), bciStartEndTimes(bciValidTrials), 'uni', 0);
binnedSpikesBci(bciValidTrials) = cellfun(@(bS) bS(:, channelsKeep), binnedSpikesBci(bciValidTrials), 'uni', 0);

binnedSpikesValidAllConcat = cat(1, binnedSpikesBci{:});
rng(0); % random number generator
numLatents = trainParams.numberFaLatents; % this is how many latents we'll project to, no matter what...
[estParams, ~] = fastfa(binnedSpikesValidAllConcat', numLatents);

latents = estParams.L; % W
[LOrthForAlignment, silde, vilde] = svd(latents);
LOrthForAlignment = LOrthForAlignment(:, 1:numLatents);


%% Check and flip the FA axis if it is needed
% doAxisFlip = false(length(trainParams.targetLatentDim),1);
doAxisFlip = false(trainParams.numberFaLatents,1);

% load unshorted electrodes from the offline dataset
% load('E:\wakko\offlineModel\offlineModel.mat')
load(['E:\wakko\offlineModel\' trainParams.offlineModelFile '.mat'])
% channelsKeepOffline = offlineDat.channelsKeep;
% loadingOffline = offlineDat.L;
channelsKeepOffline = chanKeep;
LOrthOffline = L;

% extract the common unshorted electrodes
[channelsKeepCommon, channelsKeepIndOn] = intersect(channelsKeep, channelsKeepOffline);
% channelsKeepCommonOn = channelsKeep(channelsKeepIndOn);
[~, channelsKeepIndOff] = intersect(channelsKeepOffline, channelsKeep);
% channelsKeepCommonOff = channelsKeep(channelsKeepIndOff);

% swap the FA axis based on the correlation coefficients
% axisCorrCoefs = zeros(trainParams.numberFaLatents);
% for i=1:trainParams.numberFaLatents
%     for j=1:trainParams.numberFaLatents
%         axisCorrCoef = corrcoef(latents(channelsKeepIndOn,j), loadingOffline(channelsKeepIndOff,i));
%         axisCorrCoefs(i,j) = axisCorrCoef(1,2);
%     end
% end
% 
% orthLatentsOrig = orthLatents;
% orthLatentsUpdate = orthLatents;
% estParams.Lupdate = estParams.L;
% for i=1:trainParams.numberFaLatents
%     maxCorrInd = find(abs(axisCorrCoefs(:,i))==max(abs(axisCorrCoefs(:,i))));
%     orthLatentsUpdate(:,i) = orthLatents(:,maxCorrInd);
%     estParams.Lupdate(:,i) = estParams.L(:,maxCorrInd);
% end
% estParams.L = estParams.Lupdate;
% latents = estParams.L;
% orthLatents = orthLatentsUpdate;

% flip the signs
% for i=1:length(trainParams.targetLatentDim)
if trainParams.axisSignFlip
    for i=1:trainParams.numberFaLatents
        % compute the correlation of loading using the common electrodes
        % axisCorrCoef = corrcoef(latents(channelsKeepIndOn,i), loadingOffline(channelsKeepIndOff,i));
        % if axisCorrCoef(1,2) < 0
        %     doAxisFlip(i,1) = true;
        %     estParams.L(:,i) = -estParams.L(:,i);
        %     orthLatents(:,i) = -orthLatents(:,i);
        % end

        % compute the dot product of loading using the common electrodes
        % axisCorrCoef = dot(latents(channelsKeepIndOn,i), loadingOffline(channelsKeepIndOff,i));
        axisCorrCoef = dot(LOrthForAlignment(channelsKeepIndOn,i), LOrthOffline(channelsKeepIndOff,i));
        if axisCorrCoef < 0
            doAxisFlip(i,1) = true;
            estParams.L(:,i) = -estParams.L(:,i);
            LOrthForAlignment(:,i) = -LOrthForAlignment(:,i);
        end
    end
end

if trainParams.alignAxis
    onlineLOrth = LOrthForAlignment(channelsKeepIndOn,:);
    offlineLOrth = LOrthOffline(channelsKeepIndOff,:);
    onoffLOrth = offlineLOrth' * onlineLOrth;
    [UU, DD, VV] = svd(onoffLOrth);
    O = UU * VV';
    onlineLOrthAligned = LOrthForAlignment * O';
    estParams.L = onlineLOrthAligned;
end

%% Compute the posterior mean
% Find pre uStim spike activity
% bciStimValidTrials = ~cellfun('isempty', bciEndIndex) & ~cellfun('isempty', bciStartIndex) & ~cellfun('isempty', bciStartIndex);
binnedSpikesPre = cellfun(...
        @(trialSpikeTimes, chanNums, preStartEndTime)... all the inputs
        ... adding binSizeSamples ensures that the last point is used
        histcounts2(trialSpikeTimes, chanNums, preStartEndTime(1):binSizeSamples:preStartEndTime(end)+binSizeSamples/2, 1:max(channelsKeep)+1),... histcounts2 bins things in 50ms bins and groups by channel; note that ending on the last time point means we might miss the last bin, but that's fine here
        spikeTimes(bciValidTrials), spikeChannels(bciValidTrials), preStartEndTimes(bciValidTrials), 'uni', 0);
binnedSpikesPre = cellfun(@(bS) bS(:, channelsKeep), binnedSpikesPre, 'uni', 0);
binnedSpikesPreConcat = cat(1, binnedSpikesPre{:});

% Find post uStim spike activity
binnedSpikesPost = cellfun(...
        @(trialSpikeTimes, chanNums, postStartEndTime)... all the inputs
        ... adding binSizeSamples ensures that the last point is used
        histcounts2(trialSpikeTimes, chanNums, postStartEndTime(1):binSizeSamples:postStartEndTime(end)+binSizeSamples/2, 1:max(channelsKeep)+1),... histcounts2 bins things in 50ms bins and groups by channel; note that ending on the last time point means we might miss the last bin, but that's fine here
        spikeTimes(bciValidTrials), spikeChannels(bciValidTrials), postStartEndTimes(bciValidTrials), 'uni', 0);
binnedSpikesPost = cellfun(@(bS) bS(:, channelsKeep), binnedSpikesPost, 'uni', 0);
binnedSpikesPostConcat = cat(1, binnedSpikesPost{:});

if trainParams.alignAxis
    [zPreOrth, ~, ~] = fastfa_estep(binnedSpikesPreConcat', estParams);
    [zPostOrth, ~, ~] = fastfa_estep(binnedSpikesPostConcat(1:5:end,:)', estParams);
    posteriorPre = zPreOrth.mean();
    posteriorPost = zPostOrth.mean();
    LOrth = estParams.L;
else
    % Compute posterior: 
    [ZPre, ~, beta] = fastfa_estep(binnedSpikesPreConcat', estParams);
    [zPreOrth, LOrth, TT] = orthogonalize(ZPre.mean, estParams.L);
    posteriorPre = zPreOrth; % Zdim x N trial
    %posteriorPre = ZPre.mean; % Zdim x N trial
    %zildePre = silde(1:numLatents, 1:numLatents) * vilde' * posteriorPre; % numLatents x N, notation from Byron Yu

    [ZPost, ~, beta] = fastfa_estep(binnedSpikesPostConcat(1:5:end,:)', estParams);
    [zPostOrth, LOrth, TT] = orthogonalize(ZPost.mean, estParams.L);
    posteriorPost = zPostOrth; % Zdim x N trial
    %posteriorPost = ZPost.mean; % Zdim x N trial
    %zildePost = silde(1:numLatents, 1:numLatents) * vilde' * posteriorPost; % numLatents x N, notation from Byron Yu
end

% find zilde for each uStim pattern:
stimChanCellPerTrial = cellfun(@(dS) regexp(dS, 'xippmexStimChanUpdated*=(\d*\s*\d*);*', 'tokens'), {trimmedDat.text});
% Changed ex code to send xippmexStimChanUpdated value in all types of
% trials so don't need to use the below if statement
% if trainParams.numElec > 1 % use more than 1 elec
%     stimChanCellPerTrial = cellfun(@(dS) regexp(dS, 'xippmexStimChanUpdated*=(\d*\s*\d*);*', 'tokens'), {trimmedDat.text});
% else % use only one elec
%     stimChanCellPerTrial = cellfun(@(dS) regexp(dS, 'xippmexStimChan*=(\d*);.*', 'tokens'), {trimmedDat.text});
% end

% temporary comment out the below line?
% stimChanByTrial = cellfun(@(stimChan) str2num(stimChan{1}), stimChanCellPerTrial, 'uni', 0);
stimChanByTrial = stimChanCellPerTrial;
stimChanByTrial = stimChanByTrial(bciValidTrials); % drop invalid trials
stimChanByTrial = cat(2, stimChanByTrial{:});
uniqueStimChanByTrial = unique([stimChanByTrial]);

% Maybe not necessary for closed loop codes
% stimChanByBinPre = cellfun(@(stimChan, bS) repmat(str2num(stimChan{1}), 1, size(bS, 1)), stimChanCellPerTrial(bciValidTrials), binnedSpikesPre, 'uni', 0);
% stimChanByBinPre = cat(2, stimChanByBinPre{:});
% unStimChanByBinPre = unique(stimChanByBinPre); % CHECK: unStimChanByBin or unStimChanByTrial?

% stimChanByBinPost = cellfun(@(stimChan, bS) repmat(str2num(stimChan{1}), 1, size(bS, 1)), stimChanCellPerTrial(bciValidTrials), binnedSpikesPost, 'uni', 0);
% stimChanByBinPost = cat(2, stimChanByBinPost{:});
%[unStimChanByBinPost, ~, locsUnstimChanByBinPost] = unique(stimChanByBinPost); % CHECK: unStimChanByBin or unStimChanByTrial?

targetDim = trainParams.targetLatentDim;
bciuStimChan = trainParams.bciuStimChan;
posts = zeros(length(bciuStimChan),length(targetDim));
pres = zeros(length(bciuStimChan),length(targetDim));
pushes = zeros(length(bciuStimChan),length(targetDim));
exploredCntInCalib = zeros(length(bciuStimChan),1);

% Create all uStim patterns with the number of uStim electrodes
if trainParams.numElec == 1 % singlet
   elecPatterns = cell(96, 1);
   elecPatternsInd = 1;
    for i=1:96
        elecPatterns{elecPatternsInd, 1} = num2str(i);
        elecPatternsInd = elecPatternsInd + 1;
    end
elseif trainParams.numElec == 2 % doublet
   elecPatterns = cell(4560, 1);
   % elecPatternsIndex = zeros(4560, 2);
   elecPatternsInd = 1;
   for i=1:96
       for j=1:96
           if i>=j
               continue
           end
           elecPatterns{elecPatternsInd, 1} = [num2str(i) '  ' num2str(j)];
           % elecPatterns(elecPatternsInd, 1) = i;
           % elecPatterns(elecPatternsInd, 2) = j;
           elecPatternsInd = elecPatternsInd + 1;
       end
   end
end

% Initialize the prediction with offline prediction values
% offPred = offlineDat.offPred;
postsWithOff = offPred;
exploredCntInCalibOff = zeros(length(bciuStimChan),1);

for i=1:length(uniqueStimChanByTrial)
    % if uniqueStimChanByTrial(i)==0
    if strcmp(uniqueStimChanByTrial(i), '0')||strcmp(uniqueStimChanByTrial(i), '0  0')
        continue
    end
    % curr_idx = stimChanByTrial==uniqueStimChanByTrial(i);
    curr_idx = strcmp(stimChanByTrial, uniqueStimChanByTrial(i));
    
    % update_idx = elecPatterns==uniqueStimChanByTrial(i);
    update_idx = find(matches(elecPatterns, uniqueStimChanByTrial(i)));
    curr_trials_pre = posteriorPre(targetDim, curr_idx); 
    curr_trials_post = posteriorPost(targetDim, curr_idx); 
%     posts(uniqueStimChanByTrial(i),1:length(targetDim)) = mean(curr_trials_post, 2);
%     pres(uniqueStimChanByTrial(i),1:length(targetDim)) = mean(curr_trials_pre, 2);
%     pushes(uniqueStimChanByTrial(i),1:length(targetDim)) = mean(curr_trials_post, 2) - mean(curr_trials_pre, 2);
%     exploredCntInCalib(uniqueStimChanByTrial(i),1) = sum(curr_idx);

    posts(update_idx,1:length(targetDim)) = mean(curr_trials_post, 2);
    pres(update_idx,1:length(targetDim)) = mean(curr_trials_pre, 2);
    pushes(update_idx,1:length(targetDim)) = mean(curr_trials_post, 2) - mean(curr_trials_pre, 2);
    exploredCntInCalib(update_idx,1) = sum(curr_idx);
    
    if trainParams.updateOffPred
        lr = 0.2;
        currPredError = mean(curr_trials_post, 2) - postsWithOff(update_idx,1:length(targetDim));
        postsWithOff(update_idx,1:length(targetDim)) = postsWithOff(update_idx,1:length(targetDim)) + lr*currPredError;
        
        % postsWithOff(update_idx,1:length(targetDim)) = mean(curr_trials_post, 2);
        % exploredCntInCalibOff(update_idx,1) = sum(curr_idx);
    end
end

% Initialize the ONLINE prediction table for un-tested uStim patterns
% Substitite the mean latent values for un-tested uStim patterns
% for i=targetDim
for i=1:length(targetDim)
    pres(pres(:,i)==0,i) = mean(posteriorPre(targetDim(i),:));
    posts(posts(:,i)==0,i) = mean(posteriorPost(targetDim(i),:));
    pushes(pushes(:,i)==0,i) = mean(posteriorPost(targetDim(i),:)) - mean(posteriorPre(targetDim(i),:));
end    

%%
% output: directions for task in save file
subjectCamelCase = lower(subject);
subjectCamelCase(1) = upper(subjectCamelCase(1));
bciDecoderSaveName = sprintf('%s%sFAaxisForClosedLoopStim_%s.mat', subjectCamelCase(1:2), datestr(today, 'yymmdd'), datestr(now, 'HH-MM-SS'));
%save(fullfile(bciDecoderSaveFolder, bciDecoderSaveName), 'M0', 'M1', 'M2', 'intuitiveAxisPreferredAngle', 'rotMat', 'neuralEngagementAxisInNeuralSpace', 'axisOrthToNeuralEngagementInNeuralSpace', 'channelsKeep', 'nevFilebase', 'nevFilesForTrain', 'includeBaseForTrain', 'nasNetName');
save(fullfile(bciDecoderSaveFolder, bciDecoderSaveName), 'LOrthForAlignment', 'channelsKeep', 'pres', 'posts', 'postsWithOff', 'pushes', 'exploredCntInCalib', 'exploredCntInCalibOff', 'estParams', 'binnedSpikesValidAllConcat', 'binnedSpikesPreConcat', 'binnedSpikesPostConcat', 'subject', 'LOrth', 'doAxisFlip', 'O');
decoderFileLocationAndName = fullfile(bciDecoderRelativeSaveFolder, bciDecoderSaveName);

