function newReturn = kalmanVelocityWithAnnulusDecoder(meanSpikeCount, currReturn, modelParams, expParams)

%vKF outputs
M0 = modelParams.M0;
M1 = modelParams.M1;
M2 = modelParams.M2;
if isfield(modelParams, 'rotMat')
    rotMat = modelParams.rotMat;
else
    rotMat = eye(2);
end
% currReturn is [velocityX, velocityY, initDist, initRadius, targState,
% Mapping]
currVelocity = currReturn(1:2); 
newVelocity = rotMat*(M0 + M1*rotMat'*currVelocity + M2 * meanSpikeCount);

% Internal State Axis Results
persistent currSmoothedLDAProjs

currTargetState = currReturn(5);
currAxisToUse = currReturn(6);

currMappingTargetParams = modelParams.mappingTargetParams{currAxisToUse}(currTargetState);

% Determine percentile that will be used for range based on target state
oppositePercentileForRange = expParams.oppositePercToUseForCondOne;
if currTargetState == 2
    oppositePercentileForRange= 100 - oppositePercentileForRange;
end
oppositeTargState = mod(currTargetState,2) + 1;
oppositeTargExtremePerc = prctile(modelParams.mappingTargetParams{currAxisToUse}(oppositeTargState).axisProjs, oppositePercentileForRange);

% Start of trial set initial seed value for LDA projection
if isempty(currSmoothedLDAProjs)
    % Seed initial value using percentile from opposite mapping
    currSmoothedLDAProjs = oppositeTargExtremePerc;
end
targsChangeBySD = expParams.targChangeByStd;


% initialize the range and target
if currTargetState == 1
    currTargRange = oppositeTargExtremePerc - currMappingTargetParams.mean;
    currTargReqState = currMappingTargetParams.mean - currMappingTargetParams.std*targsChangeBySD;
else
    currTargRange = currMappingTargetParams.mean - oppositeTargExtremePerc;
    currTargReqState = currMappingTargetParams.mean + currMappingTargetParams.std*targsChangeBySD;
end

% Grab decoder parameters 
multipleAxesParams = modelParams.multipleAxesParams; % num_axes x 1 struct array that contains LDA params for each axis trained in calibration
beta = modelParams.orthBeta; % Will be projection matrix to project values into FA space, 10 x neurons
estFAParams = modelParams.estFAParams;
d = estFAParams.d; % mean spike count vector during calibration trials (useful for z-scoring too)

% Zscoring parametersmodelParams.initialSeedValues(currAxisToUse, currTargetState);
zScoreSpikesMat = modelParams.zScoreSpikesMat;
zScoreSpikesMuTerm = modelParams.zScoreSpikesMuTerm;

% Exponential Smoothing parameters
alpha = expParams.alpha;

currAxisParams = multipleAxesParams(currAxisToUse);

% Z-score spikes; zscoreSpikesMat will be identity and zScoreSpikesMuTerm
% will be zero if zscoreSpikes is set to false in trainParams
zScoredSpikes = (zScoreSpikesMat * meanSpikeCount) - zScoreSpikesMuTerm;
% If FA is fitted on zScoredSpikes, d will be zero. Else, make sure to
% subtract mean.
newFaProjs = beta * (zScoredSpikes - d);
currLDAProjs = currAxisParams.projVec'*newFaProjs;

% apply exponential smoother to LDA projections
currSmoothedLDAProjs = (1-alpha)*currSmoothedLDAProjs + alpha*currLDAProjs;

% Compute how far projection is from requested internal state
% Setting small to 1 and Large to anything else
% Check if requested reward state is small else assume it is large
if currTargetState == 1
    % low
    % should be negative rewardAxisProjwhen rewardAxisProj is below smallReward Target
    newDistToTarget = currSmoothedLDAProjs - currTargReqState;
else
    % high
    % should be negative when rewardAxisProjs is above largeReward Target
    newDistToTarget = currTargReqState - currSmoothedLDAProjs;
end

% Compute distance ratio that will be used for annulus value calculation
% (r/R) in schematic
oneDimAxisRatio = newDistToTarget/currTargRange;

% Need to saturate annulus values like PFC BCI
if oneDimAxisRatio > 1
    oneDimAxisRatio = 1;
end
% If overshoots requested internal state in right direction,
% distanceToTarget will be negative
if oneDimAxisRatio < 0
    oneDimAxisRatio = 0;
end

% Set newSmoothedDist for return
newAnnulusRad = round(oneDimAxisRatio * expParams.maxAnnulusRad);
newReturn = [newVelocity(1); newVelocity(2); newDistToTarget; newAnnulusRad; currTargetState; currAxisToUse];
