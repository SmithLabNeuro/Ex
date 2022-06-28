function decoderFileLocationAndName = calibrateNeuralEngagementDecoderForBci(socketsControlComm, nevFilebase, nevFilesForTrain, trainParams, subject)

global params codes

gamma = trainParams.gamma;
channelNumbersUse = trainParams.rippleChannelNumbersInBci;
netFolder = params.nasNetFolderDataComputer;
nasNetName = trainParams.nasNetwork;
Qvalue = trainParams.kalmanQ;%100e3;
fprintf('q value is %d\n', Qvalue)

% don't double read the base if it's also a train file
if any(strcmp(nevFilesForTrain, nevFilebase))
    includeBaseForTrain = true;
else
    includeBaseForTrain = false;
end
nevFilesForTrain(strcmp(nevFilesForTrain, nevFilebase)) = [];

[nevBase,waves] = readNEV(nevFilebase);
nev = nevBase;
for nevFlInd = 1:length(nevFilesForTrain)
    [nevNx, wavesNx] = readNEV(nevFilesForTrain{nevFlInd});
    nevNx(:, 3) = nevNx(:, 3) + nev(end, 3) + 1;
    nev = [nev; nevNx];
    waves = [waves, wavesNx];
end
datBase = nev2dat(nevBase, 'nevreadflag', 1);

[slabel,nevLabelledData] = runNASNet({nev, waves},gamma, 'netFolder', netFolder, 'netname', nasNetName);
datStruct = nev2dat(nevLabelledData, 'nevreadflag', 1);

if ~includeBaseForTrain
    datStruct = datStruct(length(datBase)+1:end);
    nevLabelledData = nevLabelledData(size(nevBase, 1)+1:end, :);
end

bciDecoderSaveDrive = params.bciDecoderBasePathDataComputer;
bciDecoderRelativeSaveFolder = fullfile(subject);
bciDecoderSaveFolder = fullfile(bciDecoderSaveDrive, bciDecoderRelativeSaveFolder);
success = mkdir(bciDecoderSaveFolder);
if ~success
    fprintf('\nError creating new directory for BCI parameters\n')
    fprintf('\nkeyboard here...\n')
    keyboard
end

frThresh = trainParams.firingRateThreshold;
ffThresh = trainParams.fanoFactorThreshold;
coincThresh = trainParams.coincThresh;
coincTimeMs = trainParams.coincTimeMs;

samplingRate = params.neuralRecordingSamplingFrequencyHz; % samples/s
binSizeMs = trainParams.binSizeMs; % ms
msPerS = 1000; % 1000 ms in a second
binSizeSamples = binSizeMs/msPerS*samplingRate;

timeDecoderDelay = 0; % ms;
binDecoderDelay = timeDecoderDelay/binSizeMs;


if trainParams.orthToNeOnDecoderPlane
    try
        sendMessageWaitAck(socketsControlComm, 'kalmanDecoderForPlane');
        kalmanDecoderPath = receiveMessageSendAck(socketsControlComm);
    catch
        %     kalmanDecoderPath = 'satchel\Sa220608KalmanBci_09-55-15.mat'
        kalmanDecoderPath = 'satchel\Sa220621KalmanBci_11-55-54.mat'
    end
    kalmanDecoderFilepath = fullfile(bciDecoderSaveDrive, kalmanDecoderPath);
    modelParams = load(kalmanDecoderFilepath, 'M2', 'M1', 'M0', 'K', 'beta', 'estParams', 'channelsKeep', 'zScoreLatentMat', 'zScoreSpikesMat');
    channelsKeep = modelParams.channelsKeep;
    channelsKeepWithDig = [0 channelsKeep];
    nevLabelledData = nevLabelledData(ismember(nevLabelledData(:, 1), channelsKeepWithDig), :);
    trimmedDat = nev2dat(nevLabelledData, 'nevreadflag', 1);
else
    [trimmedDat, channelsKeep] = preprocessDat(datStruct, nevLabelledData, channelNumbersUse, binSizeMs, frThresh, ffThresh, coincThresh, coincTimeMs);
end
fprintf("\n%d channels being used\n", length(channelsKeep));

cursorPosSignal = codes.(trainParams.cursorPositionCode);
% joystickPos = cellfun(@(x) x(find(x == cursorPosSignal) + [1,2])-10000, {trimmedDat.event}, 'uni', 0);
% joystickPosTime = cellfun(@(x) x(x(:, 2) == cursorPosSignal, [3, 3]), {trimmedDat.event}, 'uni', 0);
joystickPosAndTime = cellfun(@(x)...
    [double(x(find(x == cursorPosSignal) + [1,2]))-10000,... this finds the positions (double conversion so it can go negative!)
    double(x(x(:, 2) == cursorPosSignal, 3))],... this finds the time (based on first signal)
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

if trainParams.orthToNeOnDecoderPlane
    estParams = modelParams.estParams;
%     channelsKeep = channelsKeep(ismember(channelsKeep, modelParams.channelsKeep));
%     chUsePrevLog = ismember(modelParams.channelsKeep, channelsKeep);
%     
%     estParams.minVar_bool = estParams.minVar_bool(chUsePrevLog);
%     estParams.L = estParams.L(chUsePrevLog, :);
%     estParams.Ph = estParams.Ph(chUsePrevLog);
%     estParams.d = estParams.d(chUsePrevLog);
%     modelParams.M2 = modelParams.M2(:, chUsePrevLog);
%     modelParams.beta = modelParams.beta(:, chUsePrevLog);
%     modelParams.estParams = estParams;
end

% find neural engagement axis:
binnedSpikesAll = cellfun(...
        @(trialSpikeTimes, chanNums)... all the inputs
        ... adding binSizeSamples ensures that the last point is used
        histcounts2(trialSpikeTimes, chanNums, trialSpikeTimes(1):binSizeSamples:trialSpikeTimes(end)+binSizeSamples, 1:max(channelsKeep)+1),... histcounts2 bins things in 50ms bins and groups by channel; note that ending on the last time point means we might miss the last bin, but that's fine here
        spikeTimes, spikeChannels, 'uni', 0);
binnedSpikesAll = cellfun(@(bS) bS(:, channelsKeep), binnedSpikesAll, 'uni', 0);
binnedSpikesAllConcat = cat(1, binnedSpikesAll{:});
numLatents = trainParams.numberFaLatents; % this is how many latents we'll project to, no matter what...

if ~trainParams.orthToNeOnDecoderPlane
    rng(0);
    [estParams, ~] = fastfa(binnedSpikesAllConcat', numLatents);
end

latents = estParams.L; % W
[orthLatents, silde, vilde] = svd(latents);
orthLatents = orthLatents(:, 1:numLatents);

%% training
% grab the neural activity during the BCI
bciValidTrials = ~cellfun('isempty', bciEndIndex) & ~cellfun('isempty', bciStartIndex) & trainTrials;
binnedSpikesBci = cell(size(spikeTimes));
binnedSpikesBci(bciValidTrials) = cellfun(...
        @(trialSpikeTimes, chanNums, bciStartEndTime)... all the inputs
        ... adding binSizeSamples/2 ensures that the last point is used
        histcounts2(trialSpikeTimes, chanNums, bciStartEndTime(1):binSizeSamples:bciStartEndTime(2)+binSizeSamples/2, 1:max(channelsKeep)+1),... histcounts2 bins things in 50ms bins and groups by channel
        spikeTimes(bciValidTrials), spikeChannels(bciValidTrials), bciStartEndTimes(bciValidTrials), 'uni', 0);
binnedSpikesBci(bciValidTrials) = cellfun(@(bS) bS(:, channelsKeep), binnedSpikesBci(bciValidTrials), 'uni', 0);
numBinsSmallestTrial = min(cellfun(@(x) size(x, 1), binnedSpikesBci(bciValidTrials)));  
cellOfMaxBin = num2cell(repmat(numBinsSmallestTrial, 1, length(binnedSpikesBci)));

binnedSpikesCurrStep = cell(size(spikeTimes));
binnedSpikesCurrStep(bciValidTrials) = cellfun(...
        @(bS, mxBin)...
        ... remove last binDecoderDelay from end as they don't have paired joystick data
        ... start from 2 because this is the current step (so it needs a previous!)
        bS(2:mxBin-binDecoderDelay, :),...
    binnedSpikesBci(bciValidTrials), cellOfMaxBin(bciValidTrials), 'uni', 0);
binnedSpikesPrevStep = cell(size(spikeTimes));
binnedSpikesPrevStep(bciValidTrials) = cellfun(...
        @(bS, mxBin)...
        ... remove last binDecoderDelay from end as they don't have paired joystick data
        ... also remove one more because this will be the previous step
        bS(1:mxBin-binDecoderDelay-1, :),...
    binnedSpikesBci(bciValidTrials), cellOfMaxBin(bciValidTrials), 'uni', 0);

allBinnedCountsCurrTime = cat(1, binnedSpikesCurrStep{bciValidTrials})';

allBinnedCountsCurrTimeZSc = modelParams.zScoreSpikesMat * allBinnedCountsCurrTime;
% find posterior: 
[Z, ~, beta] = fastfa_estep(allBinnedCountsCurrTimeZSc, estParams);
posterior = Z.mean; % Zdim x N

zilde = silde(1:numLatents, 1:numLatents) * vilde' * posterior; % numLatents x N, notation from Byron Yu
% zscMat = diag(1./std(zilde, [], 2));
% zilde = zscMat * zilde;

% find residuals of zilde (rezilduals):
angsCellPerTrial = cellfun(@(dS) regexp(dS, 'targetAngle*=(\d*);.*', 'tokens'), {trimmedDat.text});
angsByTrial = cellfun(@(ang) str2num(ang{1}), angsCellPerTrial, 'uni', 0);
angsByTrial = cat(2, angsByTrial);
angsByBin = cellfun(@(ang, bS) repmat(str2num(ang{1}), 1, size(bS, 1)), angsCellPerTrial, binnedSpikesCurrStep, 'uni', 0);
angsByBin = cat(2, angsByBin{:});
[unAngsByTrial, ~, locsUnAngs] = unique(angsByBin);

rezilduals = nan(size(zilde));
% angsByBin = angsForMoveSpikes(~nanTimes);
for i=1:length(unAngsByTrial)
    curr_idx = angsByBin==unAngsByTrial(i); 
    curr_trials = zilde(:, curr_idx); 
    rezilduals(:, curr_idx) = curr_trials - mean(curr_trials, 2);
    covResidualsThisAng = cov(rezilduals(:, curr_idx)', 1);
    [pcOfResidualsInFaThisAng,~,~] = svd(covResidualsThisAng);
    neuralEngagementAxisFaSpaceByAng(:, i) = pcOfResidualsInFaThisAng(:, 1);
end

% PCA
covResiduals = cov(rezilduals', 1);
[pcOfResidualsInFa,~,~] = svd(covResiduals);

% neuralEngagementAxisFaSpace = pcOfResidualsInFa(:,1); 
neuralEngagementAxisFaSpace = mean(neuralEngagementAxisFaSpaceByAng, 2);
% neuralEngagementAxisFaSpace = neuralEngagementAxisFaSpaceByAng(:, 1);
neuralEngagementAxisFaSpace = neuralEngagementAxisFaSpace/norm(neuralEngagementAxisFaSpace);

% find the other axis:
covZilde = cov(zilde', 1);
% normalize the projected axis

% We reproject the axis into high dimensional space, using the orthogonal
% basis in which it was found (remember that the posterior was reoriented
% into this basis when zilde was defined)
% this was a misunderstanding: neuralEngagementAxisInNeuralSpace = (neuralEngagementAxisFaSpace' * zscMat * silde(1:numLatents, 1:numLatents) * vilde' * beta)';
neuralEngagementAxisInNeuralSpaceNoDenoise = orthLatents * neuralEngagementAxisFaSpace;
if trainParams.orthToNeOnDecoderPlane

    % this is the denoising FA does in high dimensions, before projecting
    % into the latents
    highDimDenoise = (modelParams.estParams.L * modelParams.estParams.L' + diag(modelParams.estParams.Ph))^(-1);
    
%     [orthoNormDecodePlane, ~] = qr(modelParams.M2');
%     neuralEngagementInDecodePlane = orthoNormDecodePlane(:, 1:2) * orthoNormDecodePlane(:, 1:2)' * neuralEngagementAxisInNeuralSpace;
%     neuralEngagementInDecodePlane = neuralEngagementInDecodePlane/norm(neuralEngagementInDecodePlane);
%     space_to_search = [neuralEngagementInDecodePlane, modelParams.M2'];
%     [otherAxes,~] = qr(space_to_search);
%     axisOrthToNeuralEngagementInNeuralSpace2 = otherAxes(:, 2);
  
%     [orthLatentsOrig, svOrig, orthVOrig] = svd(modelParams.beta');
    [orthLatentsOrigLat, svOrigLat, orthVOrigLat] = svd(modelParams.estParams.L);
%     orthLatentsOrig = orthLatentsOrig(:, 1:numLatents);
    orthLatentsOrigLat = orthLatentsOrigLat(:, 1:numLatents);

    % project the neural engagement axis into the old FA space (this is for
    % the case where the FA space used to find it is not the FA space for
    % the decode plane, so trainParams.orthToNeOnDecoderPlane == false)
    neuralEngagementAxisOldFaSpace = orthLatentsOrigLat' * neuralEngagementAxisInNeuralSpaceNoDenoise;
    % normalize; we only care about direction, not scaling, in order to
    % find orthogonal
    neuralEngagementAxisOldFaSpace = neuralEngagementAxisOldFaSpace/norm(neuralEngagementAxisOldFaSpace);
    % thinking backwards about how we got the axis in the FA space helps
    % clarify this projection (and remember that in matrix math things will
    % go from right to left): first we project neural data into the FA
    % space (so multiply by modelParams.beta) to find the posterior; then
    % we projected the posterior into the orthonormalized latents (so
    % multiplying by orthVOrigLat and then svOrigLat); *then* we computed
    % the neural engagement axis from this projection, so we can now
    % project to neuralEngagementAxisOldFaSpace
    neuralEngagementAxisInNeuralSpace = (neuralEngagementAxisOldFaSpace' * diag(diag(svOrigLat)) * orthVOrigLat' * modelParams.beta * modelParams.zScoreSpikesMat)';


    % the original K comes from fitting beta * X, where X is the neural
    % data and beta = L'*(L*L' + Psi)^(-1), but here we're defining the
    % neural engagement axis in the space defined by
    % diag(diag(svOrigLat))*orthVOrigLat'*beta; in effect, to get the data
    % back into the original K space, where the Kalman filter was defined,
    % we'd have to left multiply it by diag(diag(1./svOrigLat)) and then
    % left multiply that by orthVOrigLat:
    % 1) diag(diag(1./svOrigLat)) * diag(diag(svOrigLat))*orthVOrigLat'*beta
    %       = orthVOrigLat'*beta
    % 2) orthVOrigLat * orthVOrigLat' * beta
    %       = beta
    % So correctedKSpace is the decode plane as defined from the
    % orthonormal basis in which the neural engagement axis is defined
    correctedKSpace = (modelParams.K*orthVOrigLat*diag(diag(1./svOrigLat(1:numLatents, 1:numLatents))))';
    [orthoNormFaDecodePlane, ~] = qr(correctedKSpace);

    neuralEngagementInFaDecodePlane = orthoNormFaDecodePlane(:, 1:2) * orthoNormFaDecodePlane(:, 1:2)' * neuralEngagementAxisOldFaSpace;
    
    neuralEngagementInFaDecodePlane = neuralEngagementInFaDecodePlane/norm(neuralEngagementInFaDecodePlane);
    space_to_searchK = [neuralEngagementInFaDecodePlane, correctedKSpace];
    [otherAxesK,~] = qr(space_to_searchK);
    axisOrthToNeuralEngagementInFA = otherAxesK(:, 2);%neuralEngagementInFaDecodePlane;%
    axisOrthToNeuralEngagementInNeuralSpace = (axisOrthToNeuralEngagementInFA' * svOrigLat(1:numLatents, 1:numLatents) * orthVOrigLat' * modelParams.beta * modelParams.zScoreSpikesMat)';
%     axisOrthToNeuralEngagementInNeuralSpace = axisOrthToNeuralEngagementInNeuralSpace/norm(axisOrthToNeuralEngagementInNeuralSpace);
%     neuralEngagementAxisInNeuralSpace = modelParams.beta' * neuralEngagementAxisFaSpace;
%     neuralEngagementAxisInNeuralSpace = neuralEngagementAxisInNeuralSpace/norm(neuralEngagementAxisInNeuralSpace);
else
    [pcOfFa,~,~] = svd(covZilde);
    neuralEngagementInPcPlane = pcOfFa(:, 1:2) * pcOfFa(:, 1:2)' * neuralEngagementAxisFaSpace;
    neuralEngagementInPcPlane = neuralEngagementInPcPlane/norm(neuralEngagementInPcPlane);
    space_to_search = [neuralEngagementInPcPlane, pcOfFa(:, 1:2)];
    [otherAxes,~] = qr(space_to_search);
    axisOrthToNeuralEngagement = otherAxes(:,2);
    axisOrthToNeuralEngagementInNeuralSpace = pcOfFa * axisOrthToNeuralEngagement;
end


% orient axis so most neurons are positive
if sum(sign(neuralEngagementAxisInNeuralSpace)==1)<length(neuralEngagementAxisInNeuralSpace)/2
    neuralEngagementAxisInNeuralSpace = -neuralEngagementAxisInNeuralSpace;
end

% do the same for the orthogonal axis
if sum(sign(axisOrthToNeuralEngagementInNeuralSpace)==1)<length(axisOrthToNeuralEngagementInNeuralSpace)/2
    axisOrthToNeuralEngagementInNeuralSpace = -axisOrthToNeuralEngagementInNeuralSpace;
end


%%
% allBinnedCountsCurrTime = cat(1, binnedSpikesCurrStep{bciValidTrials})';

% now grab the kinematics we'll be training on
interpJoystickPos = cell(size(spikeTimes));
interpJoystickPos(bciValidTrials) = cellfun(...
    @(jPT, bciStartEndTime)...
        ... first interpolate the x values
        [interp1(jPT(:, 3), jPT(:, 1), bciStartEndTime(1)+binSizeSamples/2:binSizeSamples:bciStartEndTime(2))',...
        ... then interpolate the y values
        interp1(jPT(:, 3), jPT(:, 2), bciStartEndTime(1)+binSizeSamples/2:binSizeSamples:bciStartEndTime(2))'],...
    joystickPosAndTime(bciValidTrials), bciStartEndTimes(bciValidTrials), 'uni', 0);
interpJoystickVel(bciValidTrials) = cellfun(...
    @(interpPos)...
    [[0 0]; diff(interpPos)/(binSizeMs/msPerS)],... start with velocity 0, then find velocity in pix/s
    interpJoystickPos(bciValidTrials), 'uni', 0);

switch trainParams.velocityToCalibrateWith
    case 'actual'
        interpJoystickKin = interpJoystickVel;
    case 'intended'
        interpJoystickKin = interpJoystickVelIntended;
    case 'intendedPositive'
        interpJoystickKin = interpJoystickVelIntendedOnlyPos;
    otherwise
        error('Training parameter velocityToCalibrateWith must either be ''actual'' or ''intended''')
end

joystickKinCurrStep = cell(size(spikeTimes));
joystickKinCurrStep(bciValidTrials) = cellfun(...
    @(iJP, mxBin)...
        ... remove binDecoderDelay bins from start as they don't have paired neural data
        ... go 2 (instead of 1) forward to pair correctly with the neural signal 'current' step
        iJP(binDecoderDelay+2:mxBin, :),...
    interpJoystickKin(bciValidTrials), cellOfMaxBin(bciValidTrials), 'uni', 0);
joystickKinPrevStep = cell(size(spikeTimes));
joystickKinPrevStep(bciValidTrials) = cellfun(...
    @(iJP, mxBin)...
        ... remove binDecoderDelay bins from start as they don't have paired neural data
        ... end one before the last step because this is the 'previous' step (so it needs a next one!)
        iJP(binDecoderDelay+1:mxBin-1, :),...
    interpJoystickKin(bciValidTrials), cellOfMaxBin(bciValidTrials), 'uni', 0);


allJoystickKinCurrTime = cat(1, joystickKinCurrStep{bciValidTrials})';
allJoystickKinPrevTime = cat(1, joystickKinPrevStep{bciValidTrials})';

nanTimes = any(isnan(allBinnedCountsCurrTime), 1)...
    | any(isnan(allJoystickKinCurrTime), 1)...
    | any(isnan(allJoystickKinPrevTime), 1);

allBinnedCountsCurrTime(:, nanTimes) = [];
allJoystickKinCurrTime(:, nanTimes) = [];
allJoystickKinPrevTime(:, nanTimes) = [];

binnedCountsInNeuralEngagement = neuralEngagementAxisInNeuralSpace' * (allBinnedCountsCurrTime - estParams.d);
binnedCountsInOrthNeuralEngagement = axisOrthToNeuralEngagementInNeuralSpace' * (allBinnedCountsCurrTime - estParams.d);
angsForMoveSpikes = cellfun(@(ang, bS) repmat(str2num(ang{1}), 1, size(bS, 1)), angsCellPerTrial(bciValidTrials), binnedSpikesCurrStep(bciValidTrials), 'uni', 0);
angsForMoveSpikes = cat(2, angsForMoveSpikes{:});
angsForMoveSpikes(nanTimes) = [];
[unAngsForMoveSpikes, ~, locAngsForMoveSpks] = unique(angsForMoveSpikes);
for unAngInd = 1:length(unAngsForMoveSpikes)
    meanOrthNEAct(unAngInd) = mean(binnedCountsInOrthNeuralEngagement( angsForMoveSpikes==unAngsByTrial(unAngInd)));
    meanNEAct(unAngInd) = mean(binnedCountsInNeuralEngagement( angsForMoveSpikes==unAngsByTrial(unAngInd)));
end
% we take the amx response, instead of the max absolute response, because we
% want the direction whose positive "thinking" leads to maximal motion
[~, prefAngInd] = max(meanOrthNEAct);
preferredAng = unAngsForMoveSpikes(prefAngInd);
fprintf('intuitive axis preferred angle is %d\n', preferredAng)


% binnedCountsInNeuralEngagement = binnedCountsInNeuralEngagement(angsForMoveSpikes==preferredAng | angsForMoveSpikes == mod(preferredAng+180, 360));
% binnedCountsInOrthNeuralEngagement = binnedCountsInOrthNeuralEngagement(angsForMoveSpikes==preferredAng | angsForMoveSpikes == mod(preferredAng+180, 360));
% allJoystickKinCurrTime = allJoystickKinCurrTime(:, angsForMoveSpikes==preferredAng | angsForMoveSpikes == mod(preferredAng+180, 360));
% allJoystickKinPrevTime = allJoystickKinPrevTime(:, angsForMoveSpikes==preferredAng | angsForMoveSpikes == mod(preferredAng+180, 360));

intuitiveAxisPreferredAngle = preferredAng;%mod(180 - preferredAng, 360);%
rotMat = [cosd(intuitiveAxisPreferredAngle) -sind(intuitiveAxisPreferredAngle); sind(intuitiveAxisPreferredAngle) cosd(intuitiveAxisPreferredAngle)];
allJoystickKinCurrTimeRotToOrthNEPrefAng = rotMat' * allJoystickKinCurrTime;
allJoystickKinPrevTimeRotToOrthNEPrefAng = rotMat' * allJoystickKinPrevTime;

regressionFeaturesOrthNeuralEngagement = [ones(size(allJoystickKinCurrTime, 2), 1), allJoystickKinPrevTimeRotToOrthNEPrefAng(1, :)', binnedCountsInOrthNeuralEngagement'];
outputParamsOrthNE = regressionFeaturesOrthNeuralEngagement \ allJoystickKinCurrTimeRotToOrthNEPrefAng(1, :)';
regressionFeaturesNeuralEngagement = [ones(size(allJoystickKinCurrTime, 2), 1), allJoystickKinPrevTimeRotToOrthNEPrefAng(2, :)', binnedCountsInNeuralEngagement'];
outputParamsNE = regressionFeaturesNeuralEngagement \ allJoystickKinCurrTimeRotToOrthNEPrefAng(2, :)' ;

% name these as with the normal Kalman BCI, so we can use the same BCI
% computer code (for now at least...)
% if trainParams.orthToNeOnDecoderPlane
%     M1 = modelParams.M1;
%     M0 = modelParams.M0;
% else
%%
if trainParams.trainVelocityKalman
    projToLat = [axisOrthToNeuralEngagementInNeuralSpace neuralEngagementAxisInNeuralSpace]';
    separateLatentDimensions = trainParams.separateLatentDimensionsInKalman;
    latentKToUse = trainParams.setAllKToValueOfLatent;
    [M0, M1, M2, A, Q, C, R, K] = trainKalmanDecoder(binnedSpikesCurrStep, bciValidTrials, joystickKinCurrStep, Qvalue, estParams, projToLat, rotMat, separateLatentDimensions, latentKToUse);
%     rotMat = eye(2);
elseif trainParams.trainNEAsState
    projToLat = modelParams.beta * modelParams.zScoreSpikesMat;
    neValues = cell(size(binnedSpikesCurrStep));
    neValues2 = cell(size(binnedSpikesCurrStep));
    neValues(bciValidTrials) = cellfun(@(bS) (bS - estParams.d') * neuralEngagementAxisInNeuralSpace, binnedSpikesCurrStep(bciValidTrials), 'uni', 0);
    neValues2(bciValidTrials) = cellfun(@(bS) neuralEngagementAxisOldFaSpace' * diag(diag(svOrigLat)) * orthVOrigLat' * projToLat * (bS' - estParams.d), binnedSpikesCurrStep(bciValidTrials), 'uni', 0);
    separateLatentDimensions = false;
    stateValues = cellfun(@(jK, nE) [jK nE], joystickKinCurrStep, neValues, 'uni', 0);
    rotMat = eye(3);
    [M0, M1, M2, A, Q, C, R, K] = trainKalmanDecoder(binnedSpikesCurrStep, bciValidTrials, stateValues, Qvalue, estParams, projToLat, rotMat, separateLatentDimensions, latentKToUse);

else
    neScal = std(binnedCountsInNeuralEngagement)/std(binnedCountsInOrthNeuralEngagement);
    if trainParams.scaleNeuralEngagementAxis
        neuralEngagementM2 = sign(outputParamsNE(3)) * abs(neScal * outputParamsOrthNE(3));
        scaleFromNe = neScal * outputParamsOrthNE(3)/outputParamsNE(3);
    else
        neuralEngagementM2 = outputParamsNE(3);
        scaleFromNe = 1;
    end
    M2 = [outputParamsOrthNE(3) 0; 0 neuralEngagementM2] * [axisOrthToNeuralEngagementInNeuralSpace'; neuralEngagementAxisInNeuralSpace'];
    
    % M1Scale = diag(svd(M2)./svd(modelParams.M2));
    M1 = [outputParamsOrthNE(2) 0; 0 outputParamsNE(2)];% M1Scale * rotMat * modelParams.M1; %
    % M0 = rotMat * modelParams.M0;
    % end
    
    % M0 explicitly set to zero-mean regression velocities;
    M0 = -mean(M1*allJoystickKinPrevTimeRotToOrthNEPrefAng + M2*allBinnedCountsCurrTime, 2); %[outputParamsOrthNE(1); outputParamsNE(1)] - M2*estParams.d + [0;-158.1555];
end
allOnesDir = ones(size(M2,2), 1);
allOnesDir = allOnesDir ./ norm(allOnesDir);
allOnesVel = rotMat * M2 * allOnesDir;
allOnesVelOrigBci = modelParams.M2 * allOnesDir;

% engagePosVel = correctedKSpace' *neuralEngagementAxisOldFaSpace;
%  correctedKSpace'*svOrigLat(1:numLatents, 1:numLatents) * orthVOrigLat'*modelParams.beta*  neuralEngagementAxisInNeuralSpace
neDir = neuralEngagementAxisInNeuralSpace;
neDir = neDir./norm(neDir);
engagePosVel = rotMat * M2 * neDir;
engagePosVelOrigBci = modelParams.M2 * neDir;
%%
% output: directions for task in save file
subjectCamelCase = lower(subject);
subjectCamelCase(1) = upper(subjectCamelCase(1));
bciDecoderSaveName = sprintf('%s%sNeuralEngagementBci_%s.mat', subjectCamelCase(1:2), datestr(today, 'yymmdd'), datestr(now, 'HH-MM-SS'));
save(fullfile(bciDecoderSaveFolder, bciDecoderSaveName), 'M0', 'M1', 'M2', 'intuitiveAxisPreferredAngle', 'rotMat', 'neuralEngagementAxisInNeuralSpace', 'axisOrthToNeuralEngagementInNeuralSpace', 'channelsKeep', 'nevFilebase', 'nevFilesForTrain', 'includeBaseForTrain', 'nasNetName');
decoderFileLocationAndName = fullfile(bciDecoderRelativeSaveFolder, bciDecoderSaveName);

%%
engagePosVelNorm = engagePosVel/norm(engagePosVel);
allOnesVelNorm = allOnesVel/norm(allOnesVel);
M0Norm = M0/norm(M0);

engagePosVelOrigBciNorm = engagePosVelOrigBci/norm(engagePosVelOrigBci);
allOnesVelOrigBciNorm = allOnesVelOrigBci/norm(allOnesVelOrigBci);
M0OrigBciNorm = modelParams.M0/norm(modelParams.M0);

updateVels = rotMat*(M2*allBinnedCountsCurrTime + M0);
figure;scatter(updateVels(1, :), updateVels(2, :), [], locAngsForMoveSpks, 'filled')
axH = gca();hold on;
xdir = [1, 0];
ydir = [0, 1];
xEng = min([abs(axH.XLim), abs(axH.YLim)]) * xdir * engagePosVelNorm;
yEng = min([abs(axH.XLim), abs(axH.YLim)]) * ydir * engagePosVelNorm;
engX = scatter(xEng, yEng, 140, 'x', 'MarkerEdgeColor',[1, 0, 0], 'LineWidth', 3);
xAO = min([abs(axH.XLim), abs(axH.YLim)]) * xdir * allOnesVelNorm;
yAO = min([abs(axH.XLim), abs(axH.YLim)]) * ydir * allOnesVelNorm;
aoX = scatter(xAO, yAO, 140, 'x', 'MarkerEdgeColor',[1, 0, 1], 'LineWidth', 3);
xM0 = min([abs(axH.XLim), abs(axH.YLim)]) * xdir * M0Norm;
yM0 = min([abs(axH.XLim), abs(axH.YLim)]) * ydir * M0Norm;
m0X = scatter(xM0, yM0, 140, 'x', 'MarkerEdgeColor',[0, 1, 1], 'LineWidth', 3);

legend([engX, aoX, m0X], 'engagement dir', 'all ones dir', 'M0 velocity');

% hold on; scatter(engagePosVel(1), engagePosVel(2), 140, 'x', 'MarkerFaceColor',[1, 0, 0]);
% hold on; scatter(allOnesVel(1), allOnesVel(2), 140, 'x', 'MarkerFaceColor',[1, 0, 0]);
axis equal
title(sprintf('velocity updates (just from spikes and offset), prefAng %d', preferredAng));

totalVels = rotMat*(M0 + M1*allJoystickKinPrevTimeRotToOrthNEPrefAng + M2*allBinnedCountsCurrTime);
figure;scatter(totalVels(1, :), totalVels(2, :), [], locAngsForMoveSpks, 'filled')
axis equal
title(sprintf('new velocities (with known prior velocities; the regression estimate, prefAng %d', preferredAng));

figure;
% binnedCountsInNeuralEngagement = neuralEngagementAxisFaSpace' * (rezilduals);
% binnedCountsInOrthNeuralEngagement = axisOrthToNeuralEngagementInFA' * (rezilduals);
angsForPlot = angsForMoveSpikes;
axONE = subplot(2,1,1);
scatter(angsForPlot, binnedCountsInOrthNeuralEngagement); hold on
plot(unAngsByTrial, meanOrthNEAct)
plot(unAngsByTrial, meanNEAct, '--')
oNEYlim = axONE.YLim;
title(sprintf('orth to neural engagement tuning (projections of binned counts during movement), prefAng %d', preferredAng))
axNE = subplot(2,1,2);
scatter(angsForPlot, binnedCountsInNeuralEngagement); hold on
plot(unAngsByTrial, meanNEAct)
plot(unAngsByTrial, meanOrthNEAct, '--')
title(sprintf('neural engagement tuning (projections of binned counts during movement), prefAng %d', preferredAng))
NEYlim = axNE.YLim;

axONE.YLim = [min([oNEYlim,NEYlim]) max([oNEYlim, NEYlim])];
axNE.YLim = [min([oNEYlim,NEYlim]) max([oNEYlim, NEYlim])];
%%
angs = cellfun(@(dS) regexp(dS, 'targetAngle*=(\d*);.*', 'tokens'), {trimmedDat.text});
angs = cellfun(@(a) str2num(a{1}), angs);

unAngs = unique(angs);
cols = jet(length(unAngs));
figure;
p(1) = subplot(2,1,1);title(sprintf('predicted trajectory, prefAng %d', preferredAng));hold on
p(2) = subplot(2,1,2);title('true trajectory');hold on;
outVelAll = {};
angsByVelAll = {};
for angInd = 1:length(unAngs)
%     trajIndBciValid = round(sum(bciValidTrials)*rand);
    allAngTraj = find(angs(bciValidTrials) == unAngs(angInd));
    for trajAngInd = 1:length(allAngTraj)
        validTrajInd = find(bciValidTrials);
        trajInd = validTrajInd(allAngTraj(trajAngInd));
        
%         trajInd = allAngTraj(trajAngInd);
        spikesTraj = (binnedSpikesCurrStep{trajInd}');
        trueTraj = joystickKinCurrStep{trajInd}';
        nanTraj = any(isnan(trueTraj), 1);
        muCurrGivenCurr = zeros(size(M2, 1), 1);%nanmean(allJoystickKinCurrTime, 2);
        outVel = zeros(size(M2, 1));%zeros(size(trueTraj(:, ~nanTraj)));
        for t = find(~nanTraj,1,'first'):size(spikesTraj, 2)
            muCurrGivenCurr = rotMat*(M1*rotMat'*muCurrGivenCurr + M2*(spikesTraj(:, t))+ M0);
            outVel(:,t) = muCurrGivenCurr;
        end
        outVelAll = [outVelAll outVel];
        angsByVelAll = [angsByVelAll unAngs(angInd)*ones(1, size(outVel, 2))];
        subplot(2,1,1);plot([cumsum(outVel(1, :)* binSizeMs/msPerS)]', [cumsum(outVel(2, :)* binSizeMs/msPerS)]', 'color', [cols(angInd, :) 1], 'LineWidth',.5);
        axis equal
        subplot(2,1,2);plot([cumsum(trueTraj(1, ~any(isnan(trueTraj), 1))* binSizeMs/msPerS)]', [cumsum( trueTraj(2, ~any(isnan(trueTraj), 1))* binSizeMs/msPerS)]', 'color', [cols(angInd, :) 1], 'LineWidth', .5);
        axis equal

    end
end

engVels = [engagePosVelNorm engagePosVelOrigBciNorm];
aoVels = [allOnesVelNorm allOnesVelOrigBciNorm];
m0Vels = [M0Norm M0OrigBciNorm];
for ind = 1:length(p)
    p1 = p(ind);
    axes(p1);
    engVel = engVels(:, ind);
    aoVel = aoVels(:, ind);
    m0Vel = m0Vels(:, ind);
    p1.XLim;
    xdir = [1, 0];
    ydir = [0, 1];
    xEng = min([abs(p1.XLim), abs(p1.YLim)]) * xdir * engVel;
    yEng = min([abs(p1.XLim), abs(p1.YLim)]) * ydir * engVel;
    engX = scatter(xEng, yEng, 140, 'x', 'MarkerEdgeColor',[1, 0, 0], 'LineWidth', 3);
    xAO = min([abs(p1.XLim), abs(p1.YLim)]) * xdir * aoVel;
    yAO = min([abs(p1.XLim), abs(p1.YLim)]) * ydir * aoVel;
    aoX = scatter(xAO, yAO, 140, 'x', 'MarkerEdgeColor',[1, 0, 1], 'LineWidth', 3);
    xM0 = min([abs(p1.XLim), abs(p1.YLim)]) * xdir * m0Vel;
    yM0 = min([abs(p1.XLim), abs(p1.YLim)]) * ydir * m0Vel;
    m0X = scatter(xM0, yM0, 140, 'x', 'MarkerEdgeColor',[0, 1, 1], 'LineWidth', 3);
    
    
    legend([engX, aoX, m0X], 'engagement dir', 'all ones dir', 'M0 velocity');
end


% keyboard