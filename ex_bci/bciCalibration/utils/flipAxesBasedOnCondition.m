function [condOneProjs, condTwoProjs, condOneTarget, condTwoTarget, condOneRange,condTwoRange, ldaParams] = flipAxesBasedOnCondition(condOneProjs, condTwoProjs,ldaParams, targChangeByStd)
    % SDs will be used to find new target
    condTwoProjsSD = std(condTwoProjs);
    condOneProjsSD = std(condOneProjs);
    % Have targets be defined
    condTwoTarget = mean(condTwoProjs);
    condOneTarget= mean(condOneProjs);
    % Flip axis only if condTwo projection is higher than condTwo projection
    if condOneTarget > condTwoTarget
        ldaParams.projVec = ldaParams.projVec*-1;
        ldaParams.projData = ldaParams.projData*-1;
        condOneProjs = -condOneProjs;
        condTwoProjs = -condTwoProjs;
    end
    % Have targets be pushed up/down by SDs
    condTwoTarget = mean(condTwoProjs) + targChangeByStd*condTwoProjsSD;
    condOneTarget= mean(condOneProjs) - targChangeByStd*condOneProjsSD;
    fprintf('Target One is %d\n', condOneTarget)
    fprintf('Target Two is %d\n', condTwoTarget)
    % Set different R values for the two conditions 
    condTwoRange = condTwoTarget - prctile(condOneProjs, 1);
    condOneRange = prctile(condTwoProjs, 99) - condOneTarget;
    fprintf('Range for cond One is %d\n', condOneRange)
    fprintf('Range for cond Two is %d\n', condTwoRange)
end