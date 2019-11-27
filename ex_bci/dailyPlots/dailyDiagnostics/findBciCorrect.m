function [ corrIdx, stepsToCorr ] = findBciCorrect( vals, lowerThresh, nTimeSteps )
% vals        - nTrials x nTimePoints
% lowerThresh - threshold which vals must be below
% nTimeSteps  - # of consecutive time steps to remain under lowerThresh
% corrIdx     - nTrials x 1
    
    meetCriterion = vals<lowerThresh;
    
    nTrials = size(vals,1);
    
    sumFilt = ones(1,nTimeSteps);
    corrIdx = false(nTrials,1);
    stepsToCorr = nan(nTrials,1);
    
    for i_trial = 1:size(vals,1)
        sumFiltDat = filter(sumFilt,1,double(meetCriterion(i_trial,:)));
        if sum(sumFiltDat==nTimeSteps)>0
            corrIdx(i_trial) = true;
            stepsToCorr(i_trial) = find(sumFiltDat==nTimeSteps,1,'first');
        end
    end
    
end

