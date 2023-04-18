function annulusRad = rewardAxisDecoder(meanSpikeCount, requestedRewardState, modelParams, expParams)
% TODO: Figure out what these mean
% meanspikeCount (already filtered by channelKeeps it seems)
% modelParams
% expParams

persistent smoothedRewardAxisProjValue

% Grab decoder parameters 
ldaParams = modelParams.ldaParams;
beta = modelParams.beta; % Will be projection matrix to project values into FA space
estFAParams = modelParams.estFAParams;
zScoreSpikesMat = modelParams.zScoreSpikesMat;
zScoreSpikesMuTerm = modelParams.zScoreSpikesMuTerm;
d = modelParams.estFAParams.d; % mean spike count vector during calibration trials (useful for z-scoring too)

% Exponential Smoothing parameters
alpha = expParams.alpha;

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
    
else
end
    % Need to compute distance 
rewardAxisDistance = 1;

% Compute smoothed distance values
newSmoothedValue = (1-alpha)*smoothedRewardAxisProjValue + alpha*rewardAxisDistance;

% Need to saturate annulus values like PFC BCi

% Set smoothedRewardAxisProjValue to newly computed value
smoothedRewardAxisProjValue = newSmoothedValue;
annulusRad = newVelocity;