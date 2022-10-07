function decoderFileLocationAndName = calibrateKalmanVelocityDecoder(~, nevFilebase, nevFilesForTrain, subject)


gamma = 0.2;
channelNumbersUse = 257:384;
netFolder = 'C:\Users\rigmdata\spikesort\nasnet\networks';
nasNetName = 'UberNet_N50_L1';

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
datBase = nev2dat(nevBase, 'nevreadflag', true);

[slabel,nevLabelledData] = runNASNet({nev, waves},gamma, 'netFolder', netFolder, 'netname', nasNetName, 'labelSpikesAsWithWrite', true);
datStruct = nev2dat(nevLabelledData, 'nevreadflag', true);

if ~includeBaseForTrain
    datStruct = datStruct(length(datBase)+1:end);
end

bciDecoderSaveDrive = params.bciDecoderBasePathDataComputer;
bciDecoderRelativeSaveFolder = fullfile(subject, datestr(today, 'yymmdd'));
bciDecoderSaveFolder = fullfile(bciDecoderSaveDrive, bciDecoderRelativeSaveFolder);
success = mkdir(bciDecoderSaveFolder);
if ~success
    fprintf('\nError creating new directory for BCI parameters\n')
    fprintf('\nkeyboard here...\n')
    keyboard
end

frThresh = 1; % Hz
ffThresh = 8; % var/mean, Hz
coincThresh = 0.2; % 20% coincident channels

samplingRate = 30000; % samples/s
binSize = 50; % ms
msPerS = 1000; % 1000 ms in a second
binSizeSamples = binSize/msPerS*samplingRate;

timeDecoderDelay = 0; % ms;
binDecoderDelay = timeDecoderDelay/binSize;


[trimmedDat, channelsKeep] = preprocessDat(datStruct, nevLabelledData, channelNumbersUse, frThresh, ffThresh, coincThresh);

% joystickPos = cellfun(@(x) x(find(x == 142) + [1,2])-10000, {trimmedDat.event}, 'uni', 0);
% joystickPosTime = cellfun(@(x) x(x(:, 2) == 142, [3, 3]), {trimmedDat.event}, 'uni', 0);
joystickPosAndTime = cellfun(@(x)...
    [double(x(find(x == 142) + [1,2]))-10000,... this finds the positions (double conversion so it can go negative!)
    double(x(x(:, 2) == 142, 3))],... this finds the time (based on first signal)
    {trimmedDat.event}, 'uni', 0);

% this is the code after which the BCI controls the cursor; could be
% FIX_OFF, or TARG_OFF, or any of the existing codes
neuralSignalForBciStart = [135];% [CURSOR_ON] expParams.bciStartCode;
% neural signal no longer controls BCI
neuralSignalForBciEnd = [136 150];% [CURSOR_OFF CORRECT] expParams.bciEndCodes;

[bciStartIndex, ~] = cellfun(@(x) find(x(:, 2)==neuralSignalForBciStart), {trimmedDat.event}, 'uni', 0);
[bciEndIndexInit, ~] = cellfun(@(x) find(x(:, 2)==neuralSignalForBciEnd), {trimmedDat.event}, 'uni', 0);
bciEndIndex = cellfun(@(x) min(x), bciEndIndexInit, 'uni', 0);
correctTrials = cellfun(@(x) any(x(:, 2)==150), {trimmedDat.event});

bciStartEndTimes = cellfun(@(evtArr, stInd, endInd) double(evtArr([stInd endInd], 3)), {trimmedDat.event}, bciStartIndex, bciEndIndex, 'uni', 0);


spikeTimes = cellfun(@(firstSpike, spikeTimeDiffs) [firstSpike; firstSpike+cumsum(uint64(spikeTimeDiffs))], {trimmedDat.firstspike}, {trimmedDat.spiketimesdiff}, 'uni', 0);
spikeChannels = cellfun(@(spikeInfo) spikeInfo(:, 1), {trimmedDat.spikeinfo}, 'uni', 0);

% this is for running FA to project the data--no need to just use the BCI
% time points, right?
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
numLatents = 10; % this is how many latents we'll project to, no matter what...
[estParams, ~] = fastfa(binnedSpikesAllConcat', numLatents);
[~, ~, beta] = fastfa_estep(binnedSpikesAllConcat', estParams);

bciValidTrials = ~cellfun('isempty', bciEndIndex) & ~cellfun('isempty', bciStartIndex) & correctTrials;
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
    [[0 0]; diff(interpPos)/(binSize/msPerS)],... start with velocity 0, then find velocity in pix/s
    interpJoystickPos(bciValidTrials), 'uni', 0);
%%
interpJoystickKin = interpJoystickVel;
joystickKinMean = nanmean(cat(1, interpJoystickKin{:}), 1)'; % remember to transpose to 2 x 1

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
allLatentProj = allLatentProj - meanLatProj;
allJoystickKinCurrTime = allJoystickKinCurrTime - joystickKinMean;
allJoystickKinPrevTime = allJoystickKinPrevTime - joystickKinMean;

% compute Kalman parameters for the state model
A = allJoystickKinCurrTime * allJoystickKinPrevTime' / (allJoystickKinPrevTime * allJoystickKinPrevTime');
Q = 1/TallLessOne*(allJoystickKinCurrTime - A*allJoystickKinPrevTime) * (allJoystickKinCurrTime - A*allJoystickKinPrevTime)';

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
M2 = K * beta;
% baseline takes care of accounting for mean offsets in training
M0 = -M1 * joystickKinMean - M2*estParams.d - K*meanLatProj;

bciDecoderSaveName = sprintf('KalmanBci_%s.mat', datestr(now, 'HH-MM-SS'));
save(fullfile(bciDecoderSaveFolder, bciDecoderSaveName), 'M0', 'M1', 'M2', 'channelsKeep', 'A', 'Q', 'C', 'R', 'nevFilebase', 'nevFilesForTrain', 'includeBaseForTrain', 'nasNetName');
decoderFileLocationAndName = fullfile(bciDecoderRelativeSaveFolder, bciDecoderSaveName);

%%
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
        muCurrGivenPrev = nanmean(allJoystickKinCurrTime, 2);
        outVel = zeros(size(trueTraj(:, ~nanTraj)));
        for t = find(~nanTraj,1,'first'):size(spikesTraj, 2)
            muCurrGivenCurr = muCurrGivenPrev + K*(spikesTraj(:, t) - C*muCurrGivenPrev);
            outVel(:,t) = muCurrGivenCurr;
            muCurrGivenPrev = A*(muCurrGivenCurr - joystickKinMean);
        end
        p1 = subplot(2,1,1);title('predicted trajectory');hold on;plot([cumsum(outVel(1, :))]', [cumsum(outVel(2, :))]', 'color', [cols(angInd, :) 1], 'LineWidth',1);
        p2 = subplot(2,1,2);title('true trajectory');hold on;plot([cumsum(trueTraj(1, ~any(isnan(trueTraj), 1)))]', [cumsum( trueTraj(2, ~any(isnan(trueTraj), 1)))]', 'color', [cols(angInd, :) 1], 'LineWidth', .5);

    end
end
% keyboard