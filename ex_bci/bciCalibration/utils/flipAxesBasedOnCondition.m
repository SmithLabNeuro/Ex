function [condOneProjs, condTwoProjs, condOneMean, condTwoMean, condOneProjsSD,condTwoProjsSD, axisParams] = flipAxesBasedOnCondition(condOneProjIndices, condTwoProjIndices, axisParams)
    % Flip axis only if condTwo projection is higher than condTwo projection
    condOneProjs = axisParams.projData(condOneProjIndices);
    condTwoProjs = axisParams.projData(condTwoProjIndices);
    if mean(condOneProjs) > mean(condTwoProjs)
        axisParams.projVec = axisParams.projVec*-1;
        axisParams.projData = axisParams.projData*-1;
        condOneProjs = -condOneProjs;
        condTwoProjs = -condTwoProjs;
    end
    % SDs will be used to find new target
    condTwoProjsSD = std(condTwoProjs);
    condOneProjsSD = std(condOneProjs);
    % Have targets be defined
    condOneMean= mean(condOneProjs);
    condTwoMean = mean(condTwoProjs);
end