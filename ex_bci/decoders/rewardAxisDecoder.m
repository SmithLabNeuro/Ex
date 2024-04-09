function newReturn = rewardAxisDecoder(meanSpikeCount, currReturn, modelParams, expParams)
% meanspikeCount (already filtered by channelKeeps): num_channels x 1
% vector
% modelParams - info saved by calibration function
% expParams - from bci_rewardAxisDecoder.xml

persistent currSmoothedOneDimAxisProjs shamReturnBinIdx shamTrialReturnIdx

currTrialRewardState = currReturn(3); % should be a reward state (ie, 1 or 3)
currRewardStructIdx = currTrialRewardState;
% set it to the right index if reward idx is 3; code below doesn't matter
% for sham trials
if currTrialRewardState == 3
    currRewardStructIdx = 2;
end

% Currently set trialTypeIdx greater than or equal to 4 to be sham trial
if currTrialRewardState >= 4
    % do sham
    if (currTrialRewardState - 3) == 1
        % select small
        potentialShamReturnVals = modelParams.shamSmallTrialReturnVals;
    else
        % select large
        potentialShamReturnVals = modelParams.shamLargeTrialReturnVals;
    end
    % Set the Bin index that you will select along with the sham trial you
    % are showing
    if isempty(shamReturnBinIdx)
        shamTrialReturnIdx = randi(length(potentialShamReturnVals));
        shamReturnBinIdx = 1;
    end
    % select random idx
    shamReturnVals = potentialShamReturnVals{shamTrialReturnIdx};
    % set return values to be the current bin of shamReturnVals
    newDispToRequestedState = shamReturnVals(1,shamReturnBinIdx);
    oneDimAxisRatio = shamReturnVals(2,shamReturnBinIdx);
    % Make sure that the indexed bin does not exceed the length of the
    % returnVals
    if shamReturnBinIdx <= (size(shamReturnVals,2) - 1)
        shamReturnBinIdx = shamReturnBinIdx + 1;
    end
else
    % Grab decoder parameters 
    orthBeta = modelParams.orthBeta; % Will be projection matrix to project values into FA space, 10 x neurons
    estFAParams = modelParams.estFAParams;
    d = estFAParams.d; % mean spike count vector during calibration trials (useful for z-scoring too)
    currRewardAxisParams = modelParams.rewardAxisParams; % key should be reward idx and value should be struct containing SD and means

    smallStateChangeBySD = expParams.smallTargChangeByStd;
    largeStateChangeBySD = expParams.largeTargChangeByStd;

    % Zscoring parameters
    zScoreSpikesMat = modelParams.zScoreSpikesMat;
    zScoreSpikesMuTerm = modelParams.zScoreSpikesMuTerm;

    % Exponential Smoothing parameters
    alpha = expParams.alpha;
    
    % meanSDStruct is a 1 x 2 struct array
    currRewardStats = currRewardAxisParams.meanSDStruct(currRewardStructIdx);

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
    if currTrialRewardState == 1
        currRequestedRewardState = currRewardStats.mean - currRewardStats.sd*smallStateChangeBySD;
    else
        % If large, get it to go above the mean
        currRequestedRewardState =  currRewardStats.mean + largeStateChangeBySD*currRewardStats.sd;
    end
    % Set the currTargRange to be 1 times the absolute value of the requested
    % state. Have initial state start from opposite end for given distribution
    initialSeedValue = currRewardStats.oppExtreme;
    % Difference from requested reward state to the opposite Extreme
    currRewardRange = abs(currRequestedRewardState - currRewardStats.oppExtreme);
    % Start of trial set initial seed value for 1D Axis projection
    if isempty(currSmoothedOneDimAxisProjs)
        currSmoothedOneDimAxisProjs = initialSeedValue; % seed the initial value to be on the other side
    end

    % Exponentially smooth LDA projection
    currSmoothedOneDimAxisProjs = (1-alpha)*currSmoothedOneDimAxisProjs + alpha*currRewardAxisProjs;

    % Determine the displacement of the current smoothed state from the requested state
    if currTrialRewardState == 1
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
end

newReturn = [newDispToRequestedState; oneDimAxisRatio];
[newReturn; currTrialRewardState]';


