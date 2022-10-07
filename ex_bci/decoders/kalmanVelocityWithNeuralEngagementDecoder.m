function newVelocityWithEngagement = kalmanVelocityWithNeuralEngagementDecoder(meanSpikeCount, currVelocityWithEngagement, modelParams, expParams)

persistent neuralEngagementValueZscoredPrev

numBinsPrevNeZscAvg = expParams.numBinsPrevNeZscAvg;
if isempty(neuralEngagementValueZscoredPrev)
    initNe = currVelocityWithEngagement(3);
    neuralEngagementValueZscoredPrev = initNe*ones(1, numBinsPrevNeZscAvg);
end

    
% Find the next velocity using a Kalman filter first
M0 = modelParams.M0;
M1 = modelParams.M1;
M2 = modelParams.M2;
if isfield(modelParams, 'rotMat')
    rotMat = modelParams.rotMat;
else
    rotMat = eye(2);
end
currVelocity = currVelocityWithEngagement(1:2);
newVelocity = rotMat*(M0 + M1*rotMat'*currVelocity + M2 * meanSpikeCount);

% Compute neural engagement value using this velocity next
interpolatedConditionMeans = modelParams.interpolatedConditionMeans;
interpolatedNeuralEngagementAxes = modelParams.interpolatedNeuralEngagementAxes;
interpolatedNeuralEngagementValueMeans = modelParams.interpolatedNeuralEngagementValueMeans;
interpolatedNeuralEngagementValueStds = modelParams.interpolatedNeuralEngagementValueStds;
d = modelParams.d;
betaOrthNorm = modelParams.betaOrthNorm;
zScoreSpikesMat = modelParams.zScoreSpikesMat;

% Find angle of velocity (from 0 to 360)
intendedAngle = atan2d(newVelocity(2),newVelocity(1));
if((intendedAngle<0))
  intendedAngle = intendedAngle + 360; % transforms domain from [-180,180] to [0,360]
end
intendedAngle = round(intendedAngle); % makes into integer
if(intendedAngle == 360)
  intendedAngle = 0;
end

% Choose the corresponding neural engagement axis
neuralEngagementAxisForCurrAng  = interpolatedNeuralEngagementAxes(:,intendedAngle+1);
conditionMeanForCurrAng         = interpolatedConditionMeans(:,intendedAngle+1);

% Project spikes into orthonormalized FA space
factorActivity = betaOrthNorm * ((zScoreSpikesMat * meanSpikeCount)  - d);

% Subtract condition-specific mean activity and project
% onto the corresponding neural engagement axis
neuralEngagementValue           = (factorActivity - conditionMeanForCurrAng)'*neuralEngagementAxisForCurrAng;
neuralEngagementValueZscored    = (neuralEngagementValue - interpolatedNeuralEngagementValueMeans(intendedAngle+1))...
  /(interpolatedNeuralEngagementValueStds(intendedAngle+1));

neuralEngagementValueZscoredTimeAvg = mean([neuralEngagementValueZscoredPrev(1:numBinsPrevNeZscAvg), neuralEngagementValueZscored]);

%                     neuralEngagementValueZscoredTimeAvg
%                     neuralEngagementValueZscoredPrev

neuralEngagementValueZscoredPrev = [neuralEngagementValueZscored neuralEngagementValueZscoredPrev(1:numBinsPrevNeZscAvg-1)];

newVelocityWithEngagement   = newVelocity;
newVelocityWithEngagement(3)   = neuralEngagementValueZscoredTimeAvg;