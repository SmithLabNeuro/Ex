function newReturn = multipleIntuitiveAxesDecoder(meanSpikeCount, currReturn, modelParams, expParams)
% meanspikeCount (already filtered by channelKeeps): num_channels x 1
% vector
% modelParams - info saved by calibration function
% expParams - from bci_rewardAxisDecoder.xml

numFaLatents = expParams.numberFaLatents;
persistent currSmoothedFAProjs

if isempty(currSmoothedFAProjs)
    currSmoothedFAProjs = zeros(numFaLatents,1);
end

currDistToTarget = currReturn(1);
currAnnulusRad = currReturn(2);
currRewardState = currReturn(3);
currAxisToUse = currReturn(4);
% Grab decoder parameters 
multipleAxesLDAParams = modelParams.multipleAxesLDAParams; % num_axes x 1 struct array that contains LDA params for each axis trained in calibration
beta = modelParams.beta; % Will be projection matrix to project values into FA space, 10 x neurons
estFAParams = modelParams.estFAParams;
d = estFAParams.d; % mean spike count vector during calibration trials (useful for z-scoring too)

% Zscoring parameters
zScoreSpikesMat = modelParams.zScoreSpikesMat;
zScoreSpikesMuTerm = modelParams.zScoreSpikesMuTerm;

% Exponential Smoothing parameters
alpha = expParams.alpha;

% Get the targets and ranges for all axes  
lowHighTargetStates = modelParams.lowHighTargetStates;
lowHighTargetRanges = modelParams.lowHighTargetRanges;

% Determine relevant parameters for mapping being currently used
condTwoTarget = modelParams.lowHighTargetStates(currAxisToUse,2);
condOneTarget = modelParams.lowHighTargetStates(currAxisToUse,1);
condTwoRange = modelParams.lowHighTargetRanges(currAxisToUse,2);
condOneRange = modelParams.lowHighTargetRanges(currAxisToUse,1);
currLDAParams = multipleAxesLDAParams(currAxisToUse);

% Z-score spikes; zscoreSpikesMat will be identity and zScoreSpikesMuTerm
% will be zero if zscoreSpikes is set to false in trainParams
zScoredSpikes = (zScoreSpikesMat * meanSpikeCount) - zScoreSpikesMuTerm;
% If FA is fitted on zScoredSpikes, d will be zero. Else, make sure to
% subtract mean.
newFaProjs = beta * (zScoredSpikes - d);

% apply exponential smoother to FA projections
newSmoothFaProj = (1-alpha) .* currSmoothedFAProjs + alpha .* newFaProjs;

% Compute LDA Projection of smoothed FA Projs
oneDimAxisProj = currLDAParams.projVec' * newSmoothFaProj; % 1 x 1?

% Compute how far projection is from requested internal state
% Setting small to 1 and Large to anything else
% Check if requested reward state is small else assume it is large
if currRewardState == 1
    % should be negative when rewardAxisProj is below smallReward Target
    newDistToTarget = oneDimAxisProj - condOneTarget;
    oneDimAxisRange = condOneRange;
else
    % should be negative when rewardAxisProjs is above largeReward Target
    newDistToTarget = condTwoTarget - oneDimAxisProj;
    oneDimAxisRange = condTwoRange;
end

% Compute distance ratio that will be used for annulus value calculation
% (r/R) in schematic
oneDimAxisRatio = newDistToTarget/oneDimAxisRange;

% Need to saturate annulus values like PFC BCI
if oneDimAxisRatio > 1
    oneDimAxisRatio = 1;
end
% If overshoots requested internal state in right direction,
% distanceToTarget will be negative
if oneDimAxisRatio < 0
    oneDimAxisRatio = 0;
end

% Set new smoothed fa projs as curr fa Projs
currSmoothedFAProjs = newSmoothFaProj;

% Set newSmoothedDist for return
newAnnulusRad = round(oneDimAxisRatio * expParams.maxAnnulusRad);
newReturn = [newDistToTarget; newAnnulusRad; currRewardState, currAxisToUse];
