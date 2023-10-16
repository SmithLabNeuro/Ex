function newReturn = rewardAxisDecoder(meanSpikeCount, currReturn, modelParams, expParams)
% meanspikeCount (already filtered by channelKeeps): num_channels x 1
% vector
% modelParams - info saved by calibration function
% expParams - from bci_rewardAxisDecoder.xml

persistent currSmoothedOneDimAxisProjs

currRewardIdx = currReturn(3); % should be a reward state (ie, 1 or 3)
% Grab decoder parameters 
orthBeta = modelParams.orthBeta; % Will be projection matrix to project values into FA space, 10 x neurons
estFAParams = modelParams.estFAParams;
d = estFAParams.d; % mean spike count vector during calibration trials (useful for z-scoring too)
currRewardAxisParams = modelParams.rewardAxisParams; % key should be reward idx and value should be struct containing SD and means

requestedStateChangeBySD = expParams.targChangeByStd;

% Zscoring parameters
zScoreSpikesMat = modelParams.zScoreSpikesMat;
zScoreSpikesMuTerm = modelParams.zScoreSpikesMuTerm;

% Exponential Smoothing parameters
alpha = expParams.alpha;
currRewardStats = currRewardAxisParams.meanSDMap(currRewardIdx);

% Z-score spikes; zscoreSpikesMat will be identity and zScoreSpikesMuTerm
% will be zero if zscoreSpikes is set to false in trainParams
zScoredSpikes = (zScoreSpikesMat * meanSpikeCount) - zScoreSpikesMuTerm;
% If FA is fitted on zScoredSpikes, d will be zero. Else, make sure to
% subtract mean. Project onto orthonormalized factors
newFaProjs = orthBeta * (zScoredSpikes - d);
% Project onto Reward Axis (LDA, don't subtract the mean)
currRewardAxisProjs = currRewardAxisParams.projVec'*newFaProjs; % scalar projection

% Determine the requested state and the appropriate initial seed value
% If small, get it to go below the mean
if currRewardIdx == 1
    currRequestedRewardState = currRewardStats.mean - currRewardStats.sd*requestedStateChangeBySD;
else
    % If large, get it to go above the mean
    currRequestedRewardState =  currRewardStats.mean + requestedStateChangeBySD*currRewardStats.sd;
end
% Set the currTargRange to be 2 times the absolute value of the requested
% state. Have initial state start from opposite end
currRewardRange = 2*abs(currRequestedRewardState);
initialSeedValue = -1*currRequestedRewardState;

% Start of trial set initial seed value for 1D Axis projection
if isempty(currSmoothedOneDimAxisProjs)
    currSmoothedOneDimAxisProjs = initialSeedValue; % seed the initial value to be on the other side
end

% Exponentially smooth LDA projection
currSmoothedOneDimAxisProjs = (1-alpha)*currSmoothedOneDimAxisProjs + alpha*currRewardAxisProjs;

% Determine the displacement of the current smoothed state from the requested state
if currRewardIdx == 1
    newDispToRequestedState = currSmoothedOneDimAxisProjs - currRequestedRewardState;
else
    newDispToRequestedState =  currRequestedRewardState - currSmoothedOneDimAxisProjs;
end

oneDimAxisRatio = newDispToRequestedState/currRewardRange;

% If overshoots requested state, saturate the feedback
if oneDimAxisRatio > 1
    oneDimAxisRatio = 1;
end
if oneDimAxisRatio < 0
    oneDimAxisRatio = 0;
end

% Set newSmoothedDist for return
newAnnulusRad = round(oneDimAxisRatio * expParams.maxAnnulusRad);
newReturn = [newDistToTarget; newAnnulusRad; currTargetState];
