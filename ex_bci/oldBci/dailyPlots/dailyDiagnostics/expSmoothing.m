function [smoothX] = expSmoothing(X,alphaVal)
%EXPSMOOTHING Summary of this function goes here
%   Detailed explanation goes here

    [xDim,nT] = size(X);
    
    smoothX = nan(xDim,nT);
    smoothX(:,1) = X(:,1);
    for ii = 2:nT
        smoothX(:,ii) = alphaVal.*X(:,ii) + (1-alphaVal).*smoothX(:,ii-1);
    end

end

