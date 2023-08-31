function decoderFileLocationAndName = calibrateKalmanVelocityAndIntuitiveAxisAnnulusDecoderForBci(~, nevFilebase, nevFilesForTrain, trainParams, subject, offlineFlag)

global params codes
% NASNet variables
gamma = trainParams.gamma;
nasNetName = trainParams.nasNetwork;
% If you ever see an empty params, make sure to run recordEX first
netFolder = params.nasNetFolderDataComputer;
% Specify channels that will be used for BCI
channelNumbersUse = trainParams.rippleChannelNumbersInBci;

%
Qvalue = trainParams.kalmanQ;%100e3;
fprintf('q value is %d\n', Qvalue)

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
[nevBase,waves] = readNEV(nevFilebase);
nev = nevBase;
for nevFlInd = 1:length(nevFilesForTrain)
    [nevNx, wavesNx] = readNEV(nevFilesForTrain{nevFlInd});
    nevNx(:, 3) = nevNx(:, 3) + nev(end, 3) + 1;
    nev = [nev; nevNx];
    waves = [waves, wavesNx];
end
datBase = nev2dat(nevBase, 'nevreadflag', true);

[slabel,nevLabelledData] = runNASNet({nev, waves},gamma, 'netFolder', netFolder, 'netname', nasNetName, 'labelSpikesAsWithWrite', true);


datStruct = nev2dat(nevLabelledData, 'nevreadflag', true);
if ~includeBaseForTrain
    nevLabelledData = nevLabelledData(size(nevBase, 1)+1:end, :);
    % adding + 2 skips the background process trial which happens on the
    % overlap
    datStruct = datStruct(length(datBase)+2:end);%nev2dat(nevLabelledData, 'nevreadflag', true);
end
%% Initialize BCI-specific parameters
bciDecoderSaveDrive = params.bciDecoderBasePathDataComputer;
bciDecoderRelativeSaveFolder = fullfile(subject);
bciDecoderSaveFolder = fullfile(bciDecoderSaveDrive, bciDecoderRelativeSaveFolder);
success = mkdir(bciDecoderSaveFolder);
if ~success
    fprintf('\nError creating new directory for BCI parameters\n')
    fprintf('\nkeyboard here...\n')
    keyboard
end
binSizeMs = trainParams.binSizeMs; % ms
msPerS = 1000; % 1000 ms in a second

% Parameters used to identify valid channels that will be used for BCI.
frThresh = trainParams.firingRateThreshold;
ffThresh = trainParams.fanoFactorThreshold;
coincThresh = trainParams.coincThresh;
coincTimeMs = trainParams.coincTimeMs;
timeDecoderDelay = 0; % ms;
binDecoderDelay = timeDecoderDelay/binSizeMs;
%% Preprocess and train decoder using calibration data
% remove all non-calibration trials, as they don't have correct spiking
% data and will mess up FA training
nonCalibrationTrialCode = codes.BACKGROUND_PROCESS_TRIAL; % Code 250 should refer to nonCalibration Trials
nonCalibrationTrialsInds = cellfun(@(x) any(ismember(x(:, 2), nonCalibrationTrialCode)), {datStruct.event});
[trimmedDat, channelsKeep] = preprocessDat(datStruct, nevLabelledData, channelNumbersUse, binSizeMs, frThresh, ffThresh, coincThresh, coincTimeMs);
trimmedDat = trimmedDat(~nonCalibrationTrialsInds);
fprintf("\n%d channels being used\n", length(channelsKeep));

% Keep track of cursor position signal
cursorPosSignal = codes.(trainParams.cursorPositionCode);
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
% Identify calibration trial indices that have matching result codes
trainTrials = cellfun(@(x) any(ismember(x(:, 2), resultCodesForTrialsToKeep)), {trimmedDat.event});

% Identify samples that correspond to the end and beginning for epoch
samplingRate = params.neuralRecordingSamplingFrequencyHz; % samples/s
binSizeSamples = binSizeMs/msPerS*samplingRate;

bciStartEndTimes = cellfun(@(evtArr, stInd, endInd) double(evtArr([stInd endInd], 3)), {trimmedDat.event}, bciStartIndex, bciEndIndex, 'uni', 0);
% Get all spike times and spike channels from dat
spikeTimes = cellfun(@(firstSpike, spikeTimeDiffs) [firstSpike; firstSpike+cumsum(uint64(spikeTimeDiffs))], {trimmedDat.firstspike}, {trimmedDat.spiketimesdiff}, 'uni', 0);
spikeChannels = cellfun(@(spikeInfo) spikeInfo(:, 1), {trimmedDat.spikeinfo}, 'uni', 0);
%% Grab Neural activity that will be used to fit the internal state axis 
binnedSpikesAll = cellfun(...
        @(trialSpikeTimes, chanNums)... all the inputs
        ... adding binSizeSamples ensures that the last point is used
        ... while subtracting 0.1 is for a super corner case where a spike overlaps the last bin end
        ... --histcounts2 includes this in the last bin (unlike what it does for every other bin);
        ... as spike times are integers the 0.1 subtraction changes no other bins
        histcounts2(trialSpikeTimes, chanNums, (trialSpikeTimes(1):binSizeSamples:trialSpikeTimes(end)+binSizeSamples)-0.1, 1:max(channelsKeep)+1),... histcounts2 bins things in 50ms bins and groups by channel; note that ending on the last time point means we might miss the last bin, but that's fine here
        spikeTimes, spikeChannels, 'uni', 0);
binnedSpikesAll = cellfun(@(bS) bS(:, channelsKeep), binnedSpikesAll, 'uni', 0);
binnedSpikesAllConcat = cat(1, binnedSpikesAll{:});
numLatents = trainParams.numberFaLatents; % this is how many latents we'll project to, no matter what...
if trainParams.zScoreSpikes
    zScoreSpikesMat = diag(1./std(binnedSpikesAllConcat, [], 1));
else
    zScoreSpikesMat = eye(size(binnedSpikesAllConcat, 2));
end
% this is for running FA to project the data--no need to just use the BCI
% time points, right?
binnedSpikesAllConcat = binnedSpikesAllConcat * zScoreSpikesMat;
[estParams, ~] = fastfa(binnedSpikesAllConcat', numLatents);
[~, ~, beta] = fastfa_estep(binnedSpikesAllConcat', estParams);

% Identify trials that will be used for training BCI (checkes that trials
% contain start and end indices
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
binnedSpikesBci(bciValidTrials) = cellfun(@(bS) bS(:, channelsKeep)*zScoreSpikesMat, binnedSpikesBci(bciValidTrials), 'uni', 0);

   
binnedSpikesCurrStep = cell(size(spikeTimes));
binnedSpikesCurrStep(bciValidTrials) = cellfun(...
        @(bS)...
        ... remove last binDecoderDelay from end as they don't have paired joystick data
        ... start from 2 because this is the current step (so it needs a previous!)
        bS(2:end-binDecoderDelay, :),...
    binnedSpikesBci(bciValidTrials), 'uni', 0);
binnedSpikesPrevStep = cell(size(spikeTimes));
binnedSpikesPrevStep(bciValidTrials) = cellfun(...
        @(bS)...
        ... remove last binDecoderDelay from end as they don't have paired joystick data
        ... also remove one more because this will be the previous step
        bS(1:end-binDecoderDelay-1, :),...
    binnedSpikesBci(bciValidTrials), 'uni', 0);


% gotta get an estimate of the position per bin, since positions aren't
% necessarily updated per bin (or even once a bin!)--currently doing a
% linear interpolation
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

paramsAll = [trimmedDat.params];
targAngles =  arrayfun(@(prm) prm.trial.targetAngle, paramsAll);
if isfield(paramsAll(1).block, 'targetDistance')
    targDists = arrayfun(@(prm) prm.block.targetDistance, paramsAll);
else
    targDists = repmat(datStruct(1).params.block.targetDistance, size(targAngles));
end
targLocs = targDists'.*[cosd(targAngles)' sind(targAngles)'];
targLocsCell = mat2cell(targLocs, ones(size(targLocs, 1), 1), size(targLocs, 2))';
interpJoystickVelIntended(bciValidTrials) = cellfun(...
    @(interpVel, interpPos, targLoc)...
    diag(interpVel * ((targLoc-interpPos)./sqrt(sum((targLoc-interpPos).^2,2)))').*((targLoc-interpPos)./sqrt(sum((targLoc-interpPos).^2,2))),... start with velocity 0, then find velocity in pix/s
    interpJoystickVel(bciValidTrials), interpJoystickPos(bciValidTrials), targLocsCell(bciValidTrials), 'uni', 0);

interpJoystickVelIntendedOnlyPos(bciValidTrials) = cellfun(...
    @(interpVel, interpPos, targLoc)...
    max(diag(interpVel * ((targLoc-interpPos)./sqrt(sum((targLoc-interpPos).^2,2)))'), 0, 'includenan').*((targLoc-interpPos)./sqrt(sum((targLoc-interpPos).^2,2))),... start with velocity 0, then find velocity in pix/s
    interpJoystickVel(bciValidTrials), interpJoystickPos(bciValidTrials), targLocsCell(bciValidTrials), 'uni', 0);

%%
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
    @(iJP)...
        ... remove binDecoderDelay bins from start as they don't have paired neural data
        ... go 2 (instead of 1) forward to pair correctly with the neural signal 'current' step
        iJP(binDecoderDelay+2:end, :),...
    interpJoystickKin(bciValidTrials), 'uni', 0);
joystickKinPrevStep = cell(size(spikeTimes));
joystickKinPrevStep(bciValidTrials) = cellfun(...
    @(iJP)...
        ... remove binDecoderDelay bins from start as they don't have paired neural data
        ... end one before the last step because this is the 'previous' step (so it needs a next one!)
        iJP(binDecoderDelay+1:end-1, :),...
    interpJoystickKin(bciValidTrials), 'uni', 0);

% all of these will have time as a column, and value (neuron/position) as
% a row--so they're variable X time
allBinnedCountsCurrTime = cat(1, binnedSpikesCurrStep{bciValidTrials})';
allJoystickKinCurrTime = cat(1, joystickKinCurrStep{bciValidTrials})';
allBinnedCountsPrevTime = cat(1, binnedSpikesPrevStep{bciValidTrials})';
allJoystickKinPrevTime = cat(1, joystickKinPrevStep{bciValidTrials})';
Tn = cellfun(@(x) size(x, 1), binnedSpikesCurrStep(bciValidTrials));
Tall = sum(Tn);
TallLessOne = sum(Tn-1);

nanTimes = any(isnan(allBinnedCountsCurrTime), 1)...
    | any(isnan(allJoystickKinCurrTime), 1)...
    | any(isnan(allBinnedCountsPrevTime), 1)...
    | any(isnan(allJoystickKinPrevTime), 1);

allBinnedCountsCurrTime(:, nanTimes) = [];
allJoystickKinCurrTime(:, nanTimes) = [];
allBinnedCountsPrevTime(:, nanTimes) = [];
allJoystickKinPrevTime(:, nanTimes) = [];
Tall = Tall - sum(nanTimes);
TallLessOne = TallLessOne - sum(nanTimes);

allLatentProj = beta * (allBinnedCountsCurrTime - estParams.d);
meanLatProj = mean(allLatentProj, 2);
stdLatProj = std(allLatentProj, [], 2);
if trainParams.zScoreLatents
    zScoreLatentMat = diag(1./stdLatProj);
else
    zScoreLatentMat = eye(length(meanLatProj));
end
allLatentProj = zScoreLatentMat*(allLatentProj - meanLatProj);
% allJoystickKinCurrTime = allJoystickKinCurrTime - joystickKinMean;
% allJoystickKinPrevTime = allJoystickKinPrevTime - joystickKinMean;

% compute Kalman parameters for the state model
A = eye(size(allJoystickKinCurrTime, 1)); %allJoystickKinCurrTime * allJoystickKinPrevTime' / (allJoystickKinPrevTime * allJoystickKinPrevTime');
Q = [Qvalue 0; 0 Qvalue]; %1/TallLessOne*(allJoystickKinCurrTime - A*allJoystickKinPrevTime) * (allJoystickKinCurrTime - A*allJoystickKinPrevTime)';

% compute Kalman parameters for the observation model
C = allLatentProj * allJoystickKinCurrTime' / (allJoystickKinCurrTime * allJoystickKinCurrTime');
R = 1/Tall * (allLatentProj - C*allJoystickKinCurrTime) * (allLatentProj - C*allJoystickKinCurrTime)';

% let's try and converge a Kalman gain?
sigCurrGivenPrev = cov(allJoystickKinCurrTime');
muCurrGivenPrev = nanmean(allJoystickKinCurrTime, 2);


for t = 1:100
    Kcurr = sigCurrGivenPrev*C' / (C*sigCurrGivenPrev*C' + R);
    sigCurrGivenCurr = sigCurrGivenPrev - Kcurr*C*sigCurrGivenPrev;
    sigCurrGivenPrev = A*sigCurrGivenCurr*A' + Q;
    Kall{t} = Kcurr;
end
K = Kall{end};

M1 = A - K*C*A;
M2 = K * zScoreLatentMat * beta * zScoreSpikesMat;
% baseline takes care of accounting for mean offsets in training
% stateModelOffset = nanmean(cat(1, interpJoystickKin{:}), 1)'; % remember to transpose to 2 x 1
stateModelOffset = [0 0]'; % for model x_t = Ax_{t-1} + stateModelOffset
% note below that estParams.d is computed post-z-scoring, which is why it
% doesn't get multiplied by zScoreSpikesMat
M0 = (eye(2) - K*C) * stateModelOffset - K * zScoreLatentMat * beta * estParams.d - K*meanLatProj;
% M0 = -M1 * joystickKinMean - M2*estParams.d - K*meanLatProj;
%% Generate plots of Predicted trajectories based on Kalman Filter
angs = cellfun(@(dS) regexp(dS, '\w*=(\d*);.*', 'tokens'), {trimmedDat.text});
angs = cellfun(@(a) str2num(a{1}), angs);

unAngs = unique(angs);
cols = jet(length(unAngs));
figure;hold on;

for angInd = 1:length(unAngs)
%     trajIndBciValid = round(sum(bciValidTrials)*rand);
    allAngTraj = find(angs(bciValidTrials) == unAngs(angInd));
    for trajAngInd = 1:length(allAngTraj)
        validTrajInd = find(bciValidTrials);
        trajInd = validTrajInd(allAngTraj(trajAngInd));
        
%         trajInd = allAngTraj(trajAngInd);
        spikesTraj = (beta*(binnedSpikesCurrStep{trajInd}' - estParams.d) - meanLatProj);
        trueTraj = joystickKinCurrStep{trajInd}';
        nanTraj = any(isnan(trueTraj), 1);
        muCurrGivenPrev = [0;0];%nanmean(allJoystickKinCurrTime, 2);
        outVel = zeros(size(trueTraj(:, ~nanTraj)));
        velPrevM = [0;0];
        outVelM = zeros(size(trueTraj(:, ~nanTraj)));
        for t = find(~nanTraj,1,'first'):size(spikesTraj, 2)
            muCurrGivenCurr = muCurrGivenPrev + K*(spikesTraj(:, t) - C*muCurrGivenPrev);
            outVel(:,t) = muCurrGivenCurr;
            muCurrGivenPrev = A*(muCurrGivenCurr - stateModelOffset);
            
            velCurrM = M0 + M1 * velPrevM + M2 * zScoreSpikesMat^(-1) * binnedSpikesCurrStep{trajInd}(t, :)';
            outVelM(:,t) = velCurrM;
            velPrevM = velCurrM;
        end
        p1 = subplot(3,1,1);title('predicted trajectory');hold on;plot([cumsum(outVel(1, :)* binSizeMs/msPerS)]', [cumsum(outVel(2, :)* binSizeMs/msPerS)]', 'color', [cols(angInd, :) 1], 'LineWidth',1);
        p2 = subplot(3,1,2);title('true trajectory');hold on;plot([cumsum(trueTraj(1, ~any(isnan(trueTraj), 1))* binSizeMs/msPerS)]', [cumsum( trueTraj(2, ~any(isnan(trueTraj), 1))* binSizeMs/msPerS)]', 'color', [cols(angInd, :) 1], 'LineWidth', .5);
        p3 = subplot(3,1,3);title('predicted trajectory with M');hold on;plot([cumsum(outVelM(1, :)* binSizeMs/msPerS)]', [cumsum(outVelM(2, :)* binSizeMs/msPerS)]', 'color', [cols(angInd, :) 1], 'LineWidth',1);

    end
end
%% Identify internal state Axis Decoder
binnedSpikesBciConcat = cat(1, binnedSpikesBci{:});
% Do not z-score counts when defining FA on BCI trial spikes
[estFAIAParams, ~] = fastfa(binnedSpikesBciConcat', numLatents);
% Beta is simply the projection matrix into the factor space
[~, ~, betaForIA] = fastfa_estep(binnedSpikesBciConcat', estFAIAParams);
% Orthonormalize your FA loadings matrix, using the economy version of SVD
[~,D,V] = svd(estFAIAParams.L, 0);
orthBetaForIA = D*V'*betaForIA;
%% Temporally Smooth Fa Projections before using LDA/PCA
alpha = trainParams.alpha;
% Find FA projections for all trials
faProjsByTrial = cellfun(@(x) orthBetaForIA*(x' - estFAIAParams.d), binnedSpikesBciConcat, 'UniformOutput', false);
% Exponentially smooth each trial's bins, include previous bins' effects.
% Seed this way for smoothing the FA latents
initialSeedValue = mean(horzcat(faProjsByTrial{:}),2); % Should be zero vector
smoothedFaProjsByTrial = cellfun(@(x) exponentialSmoother(x, alpha, nan), faProjsByTrial, 'UniformOutput', false);
validTrialTargRepeatedLabels = [];
for k = 1:length(smoothedFaProjsByTrial)
    numBinsInTrial = size(smoothedFaProjsByTrial{k},2);
    validTrialTargRepeatedLabels = [validTrialTargRepeatedLabels,  repmat(validTrialTargetLabels(k), 1,numBinsInTrial)];
end
%% Find axes based on smoothed FA Projections
axisSelectionMethod = trainParams.axisSelectionMethod;
% target axis
% 1) subsample trials to get only up and down targets
targetsToSeparate = [90, 270; 0, 180]; % num_axes x num_targets
%targetsToSeparate = [90,270];
multipleAxesParams = [];
mappingTargetParams = {};
targChangeBySD = trainParams.targChangeByStd;
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
    if strcmpi(axisSelectionMethod, 'lda') 
        axisParams = fit_LDA(trainX, trainY);
    elseif strcmpi(axisSelectionMethod, 'pca')
        % needs to have projVec
        axisParams = struct(); % 
        % Compute target PSTHs
        uniqueLabels = unique(trainY);
        targPSTHs = zeros(length(uniqueLabels), size(trainX,2));
        for i=1:length(uniqueLabels)
            % Identify Latents for Trials that have given label
            currLabel = uniqueLabels(i);
            targPSTHs(i,:) = mean(trainX(trainY == currLabel, :));
        end
        % Compute PCA on PSTHs
        axisParams.projVec = pca(targPSTHs); % returns the first component (num_latents x 1) after computing covariance matrix and then doing svd
        axisParams.mu = mean(targPSTHs); % 1 x num_latents
        axisParams.projData = (trainX - axisParams.mu)*axisParams.projVec;
    end
    firstTargIndices = find(trainY == targetsToSeparate(k,1));
    secondTargIndices = find(trainY == targetsToSeparate(k,2));
    fprintf('Current axis for separating %i and %i\n', targetsToSeparate(k,1), targetsToSeparate(k,2))
    % Flip if necessary, firstTargProjs and secondTargProjs should have
    % same order as order in projVec
    [firstTargProjs, secondTargProjs, firstTargMean, secondTargMean, firstTargProjsSD,secondTargProjsSD, axisParams] = flipAxesBasedOnCondition(firstTargIndices, secondTargIndices, axisParams);
    targetSpecificParams = struct('angle', {}, 'mean', {}, 'std', {}, 'axisProjs', {});
    % Add struct array for first target
    targetSpecificParams(end + 1) = struct('angle', targetsToSeparate(k,1), 'mean', firstTargMean, 'std',firstTargProjsSD, 'axisProjs',  firstTargProjs);
    % Add struct array for second target
    targetSpecificParams(end + 1) = struct('angle', targetsToSeparate(k,2), 'mean', secondTargMean, 'std',secondTargProjsSD, 'axisProjs',  secondTargProjs);
    % Plot projections along this 1D intuitive axis
    figure;
    hold on
    histogram(firstTargProjs, 'binWidth', 0.25, 'FaceColor', 'b' ,'DisplayName', 'Targ1')
    histogram(secondTargProjs, 'binWidth', 0.25, 'FaceColor', 'r' ,'DisplayName', 'Targ2')
    xline(firstTargMean, '--b', 'LineWidth', 2, 'DisplayName', 'Mean of Targ 1 Projs')
    xline(secondTargMean, '--r', 'LineWidth', 2, 'DisplayName', 'Mean of Targ 2 Projs')
    xline(firstTargMean - firstTargProjsSD*targChangeBySD, '--b', 'LineWidth', 3, 'DisplayName', 'Target state for 1')
    xline(secondTargMean + secondTargProjsSD*targChangeBySD, '--r', 'LineWidth', 3, 'DisplayName', 'Target state for 2')
    hold off
    title(sprintf('Projections of calibration bins along 1D intuitive where Targ1: %i Targ2: %i', targetsToSeparate(k,1), targetsToSeparate(k,2)))
    legend()
    % Identify ranges that will be used for each state 
    multipleAxesParams = [multipleAxesParams, axisParams];
    mappingTargetParams{end+1} = targetSpecificParams;
end
%% Save model parameters in decoder
subjectCamelCase = lower(subject);
subjectCamelCase(1) = upper(subjectCamelCase(1));
if offlineFlag
    % Use provided nev filename for date
    fileBackSlashIndices = find(nevFilebase == '\');
    if length(fileBackSlashIndices) > 1
        fileBackSlashIndex = fileBackSlashIndices(end);
    else
        fileBackSlashIndex = fileBackSlashIndices;
    end
    bciDecoderSaveName = sprintf('%sKalmanVelocityFilterWithIntuitiveAxesDecoderBci_%s_offline.mat',nevFilebase(fileBackSlashIndex+1:end-4) , datestr(now, 'dd-mmm-yyyy_HH-MM-SS') );
else
    % Use current date as naming convention for online
    bciDecoderSaveName = sprintf('%s%sKalmanVelocityFilterWithIntuitiveAxesDecoderBci_%s.mat', subjectCamelCase(1:2), datestr(today, 'yymmdd'), datestr(now, 'HH-MM-SS'));
end
save(fullfile(bciDecoderSaveFolder, bciDecoderSaveName), ...
    'M0', 'M1', 'M2', 'channelsKeep', 'A', 'Q', 'C', 'R', 'beta', 'K', 'zScoreSpikesMat', ...
    'zScoreLatentMat', 'estParams', 'nevFilebase', 'nevFilesForTrain', 'includeBaseForTrain', 'nasNetName',...
    'estFAIAParams', 'orthBetaForIA', 'mappingTargetParams', 'multipleAxesParams' , 'trainParams', ...
    'binnedSpikesBciConcat');
decoderFileLocationAndName = fullfile(bciDecoderRelativeSaveFolder, bciDecoderSaveName);
fprintf('decoder file saved at : %s\n', decoderFileLocationAndName)

% keyboard