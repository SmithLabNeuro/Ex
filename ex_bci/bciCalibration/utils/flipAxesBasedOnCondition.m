function [condOneProjs, condTwoProjs, condOneMean, condTwoMean, condOneProjsSD,condTwoProjsSD, ldaParams] = flipAxesBasedOnCondition(condOneProjIndices, condTwoProjIndices, ldaParams)
    % Flip axis only if condTwo projection is higher than condTwo projection
    condOneProjs = ldaParams.projData(condOneProjIndices);
    condTwoProjs = ldaParams.projData(condTwoProjIndices);
    if mean(condOneProjs) > mean(condTwoProjs)
        ldaParams.projVec = ldaParams.projVec*-1;
        ldaParams.projData = ldaParams.projData*-1;
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