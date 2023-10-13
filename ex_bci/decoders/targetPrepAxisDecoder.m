function newReturn = targetPrepAxisDecoder(meanSpikeCount, currReturn, modelParams, expParams)
% meanspikeCount (already filtered by channelKeeps): num_channels x 1
% vector
% modelParams - info saved by calibration function
% expParams - from bci_rewardAxisDecoder.xml

persistent currSmoothedOneDimAxisProjs


currTargetState = currReturn(3); % should be a target angle in degrees (e.g, 0, 45, 180, etc.)
% Grab decoder parameters 
orthBeta = modelParams.orthBeta; % Will be projection matrix to project values into FA space, 10 x neurons
estFAParams = modelParams.estFAParams;
d = estFAParams.d; % mean spike count vector during calibration trials (useful for z-scoring too)
currTargAngleMap = modelParams.targAngleMap; % key should be target angle in degrees and value should be struct 

requestedStateChangeBySD = expParams.targChangeByStd;

% Zscoring parametersmodelParams.initialSeedValues(currAxisToUse, currTargetState);
zScoreSpikesMat = modelParams.zScoreSpikesMat;
zScoreSpikesMuTerm = modelParams.zScoreSpikesMuTerm;

% Exponential Smoothing parameters
alpha = expParams.alpha;
currTargAngParams = currTargAngleMap(currTargetState);
currTargPCParams = modelParams.axisParams;

% Z-score spikes; zscoreSpikesMat will be identity and zScoreSpikesMuTerm
% will be zero if zscoreSpikes is set to false in trainParams
zScoredSpikes = (zScoreSpikesMat * meanSpikeCount) - zScoreSpikesMuTerm;
% If FA is fitted on zScoredSpikes, d will be zero. Else, make sure to
% subtract mean. Project onto orthonormalized factors
newFaProjs = orthBeta * (zScoredSpikes - d);
% Project onto Target PCs after subtracting the mean
currTargPCProjs = currTargPCParams.projVec'*(newFaProjs - currTargPCParams.mu'); % 2 x 1 projection

% Project Targ PC projection onto the target Prep axis
currTargPrepAxis = currTargAngParams.normVec;
% Requested State should include SD change
currTargPrepReqState = dot(currTargAngParams.meanTargAngProj, currTargPrepAxis) + requestedStateChangeBySD*currTargAngParams.targPrepAxisProjSD;

% Start of trial set initial seed value for 1D Axis projection
if isempty(currSmoothedOneDimAxisProjs)
    % Seed initial value as 0
    currSmoothedOneDimAxisProjs = -1*currTargPrepReqState; % seed the initial value to be on the other side
end

% Range is determined as the distance from 0 to this value;
currTargRange = 2*currTargPrepReqState;

% Norm is just to be safe
currTargPrepAxisProj = dot(currTargAngParams.meanTargAngProj, currTargPCProjs)*currTargPrepAxis/norm(currTargPrepAxis);
currTargPrepAxisDist = sqrt(sum(currTargPrepAxisProj.^2)); % should be identical to the dot product

% apply exponential smoother to 1D Axis projections (dot product onto unit
% vector which Also happens to be the distance)
currSmoothedOneDimAxisProjs = (1-alpha)*currSmoothedOneDimAxisProjs + alpha*currTargPrepAxisDist;

% Compute how far projection is from requested internal state
newDistToTarget = currTargPrepReqState - currSmoothedOneDimAxisProjs;


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
newReturn = [newDistToTarget; newAnnulusRad; currTargetState];