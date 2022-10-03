function [decoderFileLocationAndName] = calibrateTwoStepNeuralEngagementDecoderForBci(socketsControlComm, nevFilebase, nevFilesForTrain, trainParams, subject)

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
for nevFlInd = 1:length(nevFilesForTrain)
    [nevNx, wavesNx] = readNEV(nevFilesForTrain{nevFlInd});
    kr = [kr (nev(end, 3) + 1)*30000];
    nevNx(:, 3) = nevNx(:, 3) + nev(end, 3) + 1;
    nev = [nev; nevNx];
    waves = [waves, wavesNx];
end
datBase = nev2dat(nevBase, 'nevreadflag', true);

[slabel,nevLabelledData] = runNASNet({nev, waves},gamma, 'netFolder', netFolder, 'netname', nasNetName, 'labelSpikesAsWithWrite', true);


bciDecoderSaveDrive = params.bciDecoderBasePathDataComputer;
bciDecoderRelativeSaveFolder = fullfile(subject);
bciDecoderSaveFolder = fullfile(bciDecoderSaveDrive, bciDecoderRelativeSaveFolder);
% success = mkdir(bciDecoderSaveFolder);
% if ~success
%     fprintf('\nError creating new directory for BCI parameters\n')
%     fprintf('\nkeyboard here...\n')
%     keyboard
% end


samplingRate = params.neuralRecordingSamplingFrequencyHz; % samples/s
binSizeMs = trainParams.binSizeMs; % ms
msPerS = 1000; % 1000 ms in a second
binSizeSamples = binSizeMs/msPerS*samplingRate;

timeDecoderDelay = 0; % ms;
binDecoderDelay = timeDecoderDelay/binSizeMs;


% if trainParams.orthToNeOnDecoderPlane
    try
        sendMessageWaitAck(socketsControlComm, 'kalmanDecoderForPlane');
        intuitiveKalmanDecoderPath = receiveMessageSendAck(socketsControlComm);
    catch
            intuitiveKalmanDecoderPath = 'satchel\Sa220906KalmanBci_09-43-12.mat'
%             intuitiveKalmanDecoderPath = 'satchel\Sa220817KalmanBci_10-21-31.mat'
%         intuitiveKalmanDecoderPath = 'satchel\Sa220819NeuralEngagementBci_10-35-59.mat'
    end
    kalmanDecoderFilepath = fullfile(bciDecoderSaveDrive, intuitiveKalmanDecoderPath);
    modelParams = load(kalmanDecoderFilepath, 'M2', 'M1', 'M0', 'K', 'beta', 'estParams', 'channelsKeep', 'zScoreLatentMat', 'zScoreSpikesMat');
    channelsKeep = modelParams.channelsKeep;
    channelsKeepWithDig = [0 channelsKeep];
    nevLabelledData = nevLabelledData(ismember(nevLabelledData(:, 1), channelsKeepWithDig), :);
    trimmedDat = nev2dat(nevLabelledData, 'nevreadflag', true);
% else
%     [trimmedDat, channelsKeep] = preprocessDat(datStruct, nevLabelledData, channelNumbersUse, binSizeMs, frThresh, ffThresh, coincThresh, coincTimeMs);
% end
fprintf("\n%d channels being used\n", length(channelsKeep));

if ~includeBaseForTrain
    trimmedDat = trimmedDat(length(datBase)+1:end);
    nevLabelledData = nevLabelledData(size(nevBase, 1)+1:end, :);
end

cursorPosSignal = codes.(trainParams.cursorPositionCode);
joystickPosAndTime = cellfun(@(nevEvents)...
    [double(nevEvents(size(nevEvents, 1) + find(nevEvents(:, 2) == cursorPosSignal) + [1,2]))-10000,... this finds the positions (double conversion so it can go negative!)
    double(nevEvents(nevEvents(:, 2) == cursorPosSignal, 3))],... this finds the time (based on first signal)
    {trimmedDat.event}, 'uni', 0);

% this is the code after which the BCI controls the cursor; could be
% FIX_OFF, or TARG_OFF, or any of the existing codes
neuralSignalForBciStart = codes.(trainParams.bciStartCode);% [CURSOR_ON] expParams.bciStartCode;
% neural signal no longer controls BCI
neuralSignalForBciEnd = codes.(trainParams.bciEndCode);% [CURSOR_OFF CORRECT] expParams.bciEndCodes;

[bciStartIndex, ~] = cellfun(@(x) find(x(:, 2)==neuralSignalForBciStart), {trimmedDat.event}, 'uni', 0);
[bciEndIndexInit, ~] = cellfun(@(x) find(x(:, 2)==neuralSignalForBciEnd), {trimmedDat.event}, 'uni', 0);
bciEndIndex = cellfun(@(x) min(x), bciEndIndexInit, 'uni', 0);
if iscell(trainParams.trainingResultCodes)
    resultCodesForTrialsToKeep = cellfun(@(resCode) codes.(resCode),  trainParams.trainingResultCodes);
elseif ischar(trainParams.trainingResultCodes)
    resultCodesForTrialsToKeep = codes.(trainParams.trainingResultCodes);
end
trainTrials = cellfun(@(x) any(ismember(x(:, 2), resultCodesForTrialsToKeep)), {trimmedDat.event});

bciStartEndTimes = cellfun(@(evtArr, stInd, endInd) double(evtArr([stInd endInd], 3)), {trimmedDat.event}, bciStartIndex, bciEndIndex, 'uni', 0);


spikeTimes = cellfun(@(firstSpike, spikeTimeDiffs) [firstSpike; firstSpike+cumsum(uint64(spikeTimeDiffs))], {trimmedDat.firstspike}, {trimmedDat.spiketimesdiff}, 'uni', 0);
spikeChannels = cellfun(@(spikeInfo) spikeInfo(:, 1), {trimmedDat.spikeinfo}, 'uni', 0);


% find neural engagement axis:
% binnedSpikesAll = cellfun(...
%         @(trialSpikeTimes, chanNums)... all the inputs
%         ... adding binSizeSamples ensures that the last point is used
%         ... while subtracting 0.1 is for a super corner case where a spike overlaps the last bin end
%         ... --histcounts2 includes this in the last bin (unlike what it does for every other bin);
%         ... as spike times are integers the 0.1 subtraction changes no other bins
%         histcounts2(trialSpikeTimes, chanNums, (trialSpikeTimes(1):binSizeSamples:trialSpikeTimes(end)+binSizeSamples)-0.1, 1:max(channelsKeep)+1),... histcounts2 bins things in 50ms bins and groups by channel; note that ending on the last time point means we might miss the last bin, but that's fine here
%         spikeTimes, spikeChannels, 'uni', 0);
% binnedSpikesAll = cellfun(@(bS) bS(:, channelsKeep), binnedSpikesAll, 'uni', 0);
% binnedSpikesAllConcat = cat(1, binnedSpikesAll{:});
% numLatents = trainParams.numberFaLatents; % this is how many latents we'll project to, no matter what...
% 

[orthLatentsOrigLat, svOrigLat, orthVOrigLat] = svd(modelParams.estParams.L, 'econ');

%% training
% grab the neural activity during the BCI
bciValidTrials = ~cellfun('isempty', bciEndIndex) & ~cellfun('isempty', bciStartIndex) & trainTrials;
binnedSpikesBci = cell(size(spikeTimes));
binnedSpikesBci(bciValidTrials) = cellfun(...
        @(trialSpikeTimes, chanNums, bciStartEndTime)... all the inputs
        ... adding binSizeSamples/2 ensures that the last point is used
        ... while subtracting 0.1 is for a super corner case where a spike overlaps the last bin end
        ... --histcounts2 includes this in the last bin (unlike what it does for every other bin);
        ... as spike times are integers the 0.1 subtraction changes no other bins
        histcounts2(trialSpikeTimes, chanNums, (bciStartEndTime(1):binSizeSamples:bciStartEndTime(2)+binSizeSamples/2)-0.1, 1:max(channelsKeep)+1),... histcounts2 bins things in 50ms bins and groups by channel
        spikeTimes(bciValidTrials), spikeChannels(bciValidTrials), bciStartEndTimes(bciValidTrials), 'uni', 0);
binnedSpikesBci(bciValidTrials) = cellfun(@(bS) bS(:, channelsKeep), binnedSpikesBci(bciValidTrials), 'uni', 0);

if trainParams.useSameNumberOfCalibrationBinsPerTarget
    numBinsSmallestTrial = min(cellfun(@(x) size(x, 1), binnedSpikesBci(bciValidTrials)));
    cellOfMaxBin = num2cell(repmat(numBinsSmallestTrial, 1, length(binnedSpikesBci)));
else
    cellOfMaxBin = cellfun(@(x) size(x, 1), binnedSpikesBci, 'uni', 0);
end

binnedSpikesCurrStep = cell(size(spikeTimes));
binnedSpikesCurrStep(bciValidTrials) = cellfun(...
        @(bS, mxBin)...
        ... remove last binDecoderDelay from end as they don't have paired joystick data
        ... start from 2 because this is the current step (so it needs a previous!)
        bS(2:mxBin-binDecoderDelay, :),...
    binnedSpikesBci(bciValidTrials), cellOfMaxBin(bciValidTrials), 'uni', 0);

allBinnedCountsBciValidTrials = cat(1, binnedSpikesCurrStep{bciValidTrials})';
%% now grab the kinematics we'll be training on and align to neural data
% interpJoystickPos = cell(size(spikeTimes));
% interpJoystickPos(bciValidTrials) = cellfun(...
%     @(jPT, bciStartEndTime)...
%         ... first interpolate the x values
%         [interp1(jPT(:, 3), jPT(:, 1), bciStartEndTime(1)+binSizeSamples/2:binSizeSamples:bciStartEndTime(2))',...
%         ... then interpolate the y values
%         interp1(jPT(:, 3), jPT(:, 2), bciStartEndTime(1)+binSizeSamples/2:binSizeSamples:bciStartEndTime(2))'],...
%     joystickPosAndTime(bciValidTrials), bciStartEndTimes(bciValidTrials), 'uni', 0);
% interpJoystickVel(bciValidTrials) = cellfun(...
%     @(interpPos)...
%     [[0 0]; diff(interpPos)/(binSizeMs/msPerS)],... start with velocity 0, then find velocity in pix/s
%     interpJoystickPos(bciValidTrials), 'uni', 0);
% 
% switch trainParams.velocityToCalibrateWith
%     case 'actual'
%         interpJoystickKin = interpJoystickVel;
%     otherwise
%         error('Training parameter velocityToCalibrateWith must either be ''actual'' or ''intended''')
% end
% 
% joystickKinCurrStep = cell(size(spikeTimes));
% joystickKinCurrStep(bciValidTrials) = cellfun(...
%     @(iJP, mxBin)...
%         ... remove binDecoderDelay bins from start as they don't have paired neural data
%         ... go 2 (instead of 1) forward to pair correctly with the neural signal 'current' step
%         iJP(binDecoderDelay+2:mxBin, :),...
%     interpJoystickKin(bciValidTrials), cellOfMaxBin(bciValidTrials), 'uni', 0);
% joystickKinPrevStep = cell(size(spikeTimes));
% joystickKinPrevStep(bciValidTrials) = cellfun(...
%     @(iJP, mxBin)...
%         ... remove binDecoderDelay bins from start as they don't have paired neural data
%         ... end one before the last step because this is the 'previous' step (so it needs a next one!)
%         iJP(binDecoderDelay+1:mxBin-1, :),...
%     interpJoystickKin(bciValidTrials), cellOfMaxBin(bciValidTrials), 'uni', 0);
% 
% 
% allJoystickKinCurrTime = cat(1, joystickKinCurrStep{bciValidTrials})';
% allJoystickKinPrevTime = cat(1, joystickKinPrevStep{bciValidTrials})';
% allBinnedCountsCurrTime = allBinnedCountsBciValidTrials;
% 
% nanTimes = any(isnan(allBinnedCountsCurrTime), 1)...
%     | any(isnan(allJoystickKinCurrTime), 1)...
%     | any(isnan(allJoystickKinPrevTime), 1);
% 
% allBinnedCountsCurrTime(:, nanTimes) = [];
% allJoystickKinCurrTime(:, nanTimes) = [];
% allJoystickKinPrevTime(:, nanTimes) = [];

%% find neural engagement axis - step 0,1
% modelParams.zScoreSpikesMat is the inverse std during calibration.
% modelParams.estParams.d is the (mean/std) during calibration
allBinnedCountsCurrTimeZSc = modelParams.zScoreSpikesMat * allBinnedCountsBciValidTrials; 
% find posterior: 
[Z, ~, beta] = fastfa_estep(allBinnedCountsCurrTimeZSc, modelParams.estParams);
posterior = Z.mean; % Zdim x N
zilde = svOrigLat * orthVOrigLat' * posterior; % projections into orthonormalized FA space

% find residuals of zilde (rezilduals):
angsCellPerTrial = cellfun(@(dS) regexp(dS, 'targetAngle*=(\d*);.*', 'tokens'), {trimmedDat.text});
angsByTrial = cellfun(@(ang) str2num(ang{1}), angsCellPerTrial, 'uni', 0);
angsByTrial = cat(2, angsByTrial);
[unAngsByTrial, ~, locsUnAngsByTrial] = unique([angsByTrial{:}]);

angsByBin = cellfun(@(ang, bS) repmat(str2num(ang{1}), 1, size(bS, 1)), angsCellPerTrial, binnedSpikesCurrStep, 'uni', 0);
angsByBin = cat(2, angsByBin{:});
[unAngsByTrial, ~, locsUnAngsByBin] = unique(angsByBin);

rezilduals = nan(size(zilde));
neuralEngagementValue = nan(size(angsByBin));
% angsByBin = angsForMoveSpikes(~nanTimes);
for i=1:length(unAngsByTrial)
    curr_idx = angsByBin==unAngsByTrial(i); 
    curr_trials = zilde(:, curr_idx); 
    faConditionMean(:, i) = mean(curr_trials, 2);
    rezilduals(:, curr_idx) = curr_trials - faConditionMean(:, i);
    covResidualsThisAng = cov(rezilduals(:, curr_idx)', 1);
    [pcOfResidualsInFaThisAng,~,~] = svd(covResidualsThisAng);
    neuralEngagementAxisFaSpaceByAng(:, i) = pcOfResidualsInFaThisAng(:, 1);
    
    % Flip axes to have more positive weights (in neural space)
    neuralEngagementAxisNeuralSpaceByAng = orthLatentsOrigLat * neuralEngagementAxisFaSpaceByAng(:, i);

    propPosUnits(i) = mean(sign(neuralEngagementAxisNeuralSpaceByAng));
    % NOTE Change on 7 Sept 2022 after experiment: this was <0.5, now it's
    % <0
    if(propPosUnits(i)<0)
        neuralEngagementAxisFaSpaceByAng(:,i) = -neuralEngagementAxisFaSpaceByAng(:,i);
    end
    
    % Project FA activity onto corresponding neural engagement axis for
    % each target
    neuralEngagementValue(curr_idx) = rezilduals(:, curr_idx)'*neuralEngagementAxisFaSpaceByAng(:,i);
    meanNeuralEngagementValue(i) = mean(neuralEngagementValue(curr_idx));
    stdNeuralEngagementValue(i)  = std(neuralEngagementValue(curr_idx));
end
propPosUnits
%% step 2,3 interpolation

% neuralEngagementAxisFaSpaceByAng/neFaMean are 10x8, interpolation to get 10x360

% NE axis interpolation
interpAngs = 0:359;
% ensures angles are between 0 and 360, so that we can grab the center of
% the interpolation and it's smooth
unAngsByTrial = mod(unAngsByTrial, 360);
% repeated three times because we want to take the 0-360 degree range and
% make sure it's unbroken at the edges
% (i.e. if we only interpolated 0 to 720, and then we grabbed the 0 to 360
% range, nothing ensures that the value at 0 smoothly becomes the value at
% 360; if instead we interpolate -360 to 720 and then grab the 0 to 360
% range, the value at 0 smoothly interpolates to the value at 360, because
% the value at 360 needed to interpolate 359-361, which was equivalent to
% -1 to 1 that the value at 0 interpolated)
unAngsByTrialLooped = [unAngsByTrial-360 unAngsByTrial unAngsByTrial+360];
neuralEngagementAxisFaSpaceByAngLooped = [neuralEngagementAxisFaSpaceByAng, neuralEngagementAxisFaSpaceByAng, neuralEngagementAxisFaSpaceByAng]';

nd = size(neuralEngagementAxisFaSpaceByAngLooped,2);

neuralEngagementAxisFaSpaceByAngLooped = [zeros(1,nd); neuralEngagementAxisFaSpaceByAngLooped ; zeros(1,nd)]';
neuralEngagementAxisFaSpaceInterp = spline(deg2rad(unAngsByTrialLooped),neuralEngagementAxisFaSpaceByAngLooped,deg2rad(interpAngs));

% Condition mean interpolation
neFaMeanLooped = [faConditionMean, faConditionMean, faConditionMean]';
neFaMeanLooped = [zeros(1,nd); neFaMeanLooped ; zeros(1,nd)]';
neFaMeanInterp = spline(deg2rad(unAngsByTrialLooped),neFaMeanLooped,deg2rad(interpAngs));

% Interpolating neural engagement means and stds
meanNeuralEngagementValueLooped = [meanNeuralEngagementValue, meanNeuralEngagementValue, meanNeuralEngagementValue]';
meanNeuralEngagementValueLooped = [0; meanNeuralEngagementValueLooped ; 0]';
meanNeuralEngagementValueInterp = spline(deg2rad(unAngsByTrialLooped),meanNeuralEngagementValueLooped,deg2rad(interpAngs));

stdNeuralEngagementValueLooped = [stdNeuralEngagementValue, stdNeuralEngagementValue, stdNeuralEngagementValue]';
stdNeuralEngagementValueLooped = [0; stdNeuralEngagementValueLooped ; 0]';
stdNeuralEngagementValueInterp = spline(deg2rad(unAngsByTrialLooped),stdNeuralEngagementValueLooped,deg2rad(interpAngs));

%% 

M0 = modelParams.M0;
M1 = modelParams.M1;
M2 = modelParams.M2;
d  = modelParams.estParams.d;
betaOrthNorm = svOrigLat * orthVOrigLat' * modelParams.beta;
zScoreSpikesMat = modelParams.zScoreSpikesMat;

interpolatedConditionMeans = neFaMeanInterp;
interpolatedNeuralEngagementAxes = neuralEngagementAxisFaSpaceInterp;
interpolatedNeuralEngagementValueMeans = meanNeuralEngagementValueInterp;
interpolatedNeuralEngagementValueStds  = stdNeuralEngagementValueInterp;
%%
% output: directions for task in save file
subjectCamelCase = lower(subject);
subjectCamelCase(1) = upper(subjectCamelCase(1));
bciDecoderSaveName = sprintf('%s%sNeuralEngagementBci_%s.mat', subjectCamelCase(1:2), datestr(today, 'yymmdd'), datestr(now, 'HH-MM-SS'));
save(fullfile(bciDecoderSaveFolder, bciDecoderSaveName), 'M0', 'M1', 'M2','interpolatedConditionMeans', 'interpolatedNeuralEngagementAxes' ,'interpolatedNeuralEngagementValueMeans','interpolatedNeuralEngagementValueStds', 'zScoreSpikesMat','d','betaOrthNorm', 'channelsKeep', 'nevFilebase', 'nevFilesForTrain', 'includeBaseForTrain', 'nasNetName', 'intuitiveKalmanDecoderPath');
decoderFileLocationAndName = fullfile(bciDecoderRelativeSaveFolder, bciDecoderSaveName);