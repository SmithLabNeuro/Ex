function [ tuneCurves, tuneCurveSE, condLabels ] = getTuningCurves( X, conds )
% Inputs
%    -     X: neural activity/spike counts (nNeurons x nTrials)
%    - conds: condition labels (nTrials x 1)
%
% Outputs
%    - tuneCurves: avg neuronal response to each condition (nNeurons x nConditions)
%    - tuneCurveSE: SE for each condition (nNeurons x nConditions)
%    - condLabels: ordering of condition labels for tuneCurves (nConditions x 1)
%
% Akash Umakantha (aumakant@andrew.cmu.edu)
%

    [nNeurons,~] = size(X);
    condLabels = unique(conds);
    nConditions = length(condLabels);
    
    tuneCurves = nan(nNeurons,nConditions);
    tuneCurveSE = nan(nNeurons,nConditions);
    for i_cond = 1:nConditions
        currDat = X(:,conds==condLabels(i_cond));
        tuneCurves(:,i_cond) = mean(currDat,2);
        tuneCurveSE(:,i_cond) = std(currDat,1,2)./sqrt(size(currDat,2));
    end

end

