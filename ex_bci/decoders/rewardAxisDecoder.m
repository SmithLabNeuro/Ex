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
d = modelParams.estFAParams.d; % mean spike count vector during calibration trials (useful for z-scoring too)

% Exponential Smoothing parameters
alpha = expParams.alpha;

% Project zscored spikes into  FA space
% TODO: Should we subtract the mean spike count vector
currFAProjs = beta * ((zScoreSpikesMat * meanSpikeCount)  - d);
% is d = mean or mean/std? if it's just the mean then: 
% currFAProjs = beta * ((zScoreSpikesMat * (meanSpikeCount  - d));

% Compute LDA Projection of FA Projs
rewardAxisProj = currFAProjs*ldaParams.projVec;

% Compute how far projection is from requested internal state
% Need to compute distance 
rewardAxisDistance = 1;

% Compute smoothed distance values
newSmoothedValue = (1-alpha)*smoothedRewardAxisProjValue + alpha*rewardAxisDistance;

% Need to saturate annulus values like PFC BCi

% Set smoothedRewardAxisProjValue to newly computed value
smoothedRewardAxisProjValue = newSmoothedValue;
annulusRad = newVelocity;