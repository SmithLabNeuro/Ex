function [distToTarget, annulusRad] = rewardAxisDecoder(meanSpikeCount, requestedRewardState, modelParams, expParams)
% TODO: Figure out what these mean
% meanspikeCount (already filtered by channelKeeps it seems)
% modelParams - info saved by calibration function
% expParams - from bci_rewardAxisDecoder.xml

persistent smoothedRewardAxisDist

% Grab decoder parameters 
ldaParams = modelParams.ldaParams;
beta = modelParams.beta; % Will be projection matrix to project values into FA space
estFAParams = modelParams.estFAParams;
d = estFAParams.d; % mean spike count vector during calibration trials (useful for z-scoring too)

%Zscoring parameters
zScoreSpikesMat = modelParams.zScoreSpikesMat;
zScoreSpikesMuTerm = modelParams.zScoreSpikesMuTerm;

% Exponential Smoothing parameters
alpha = expParams.alpha;

% Distance Parameters
largeRewardTarget = modelParams.largeRewardMeanProj;
smallRewardTarget = modelParams.smallRewardMeanProj;
rewardAxisRange = modelParams.rewardAxisRange;

% Z-score spikes; zscoreSpikesMat will be identity and zScoreSpikesMuTerm
% will be zero if zscoreSpikes is set to false in trainParams
zScoredSpikes = (zScoreSpikesMat * meanSpikeCount) - zScoreSpikesMuTerm;
% If FA is fitted on zScoredSpikes, d will be zero. Else, make sure to
% subtract mean.
currFAProjs = beta * (zScoredSpikes - d);
% Compute LDA Projection of FA Projs
rewardAxisProj = currFAProjs*ldaParams.projVec;

% Compute how far projection is from requested internal state
% Setting small to 1 and Large to anything else
% Check if requested reward state is small else assume it is large
if requestedRewardState == 1
    % should be negative when rewardAxisProj is below smallReward Target
    distToTarget = rewardAxisProj - smallRewardTarget;
else
    % should be negative when rewardAxisProjs is above largeReward Target
    distToTarget = largeRewardTarget - rewardAxisProj;
end

% Compute smoothed distance values
if isempty(smoothedRewardAxisDist)
    % first smoothed value
    newSmoothedDist = distToTarget;
else
    newSmoothedDist = (1-alpha)*smoothedRewardAxisDist + alpha*distToTarget;
end

% If overshoots requested internal state in right direction,
% distanceToTarget will be negative
if newSmoothedDist < 0
    newSmoothedDist = 0;
end

% Compute distance ratio that will be used for annulus value calculation
% (r/R) in schematic
rewardAxisRatio = newSmoothedDist/rewardAxisRange;

% Need to saturate annulus values like PFC BCi
if rewardAxisRatio > 1
    rewardAxisRatio = 1;
end

% Set smoothedRewardAxisProjValue to newly computed value
smoothedRewardAxisDist = newSmoothedDist;

% calculate radius from ratio r/R = a/A
annulusRad = round(rewardAxisRatio * expParams.maxAnnulusRad);