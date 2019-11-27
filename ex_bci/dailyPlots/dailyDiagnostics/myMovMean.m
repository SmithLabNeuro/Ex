function [ xSmooth ] = myMovMean( x, windowSize )
%SEPBCIDAT Summary of this function goes here
%   Detailed explanation goes here

    if ~isvector(x)
        error('Input x must be a vector');
    end

    N = length(x);
    kb = windowSize(1);
    kf = windowSize(2);
    xSmooth = nan(size(x));
    
    % smooth vector x using indexes ii-kb to ii+kf
    for ii = 1:N
        startIdx = max(ii-kb,1); % startIdx cannot be smaller than 1
        endIdx = min(ii+kf,N); % endIdx cannot be greater than N
        xSmooth(ii) = nanmean(x(startIdx:endIdx));
    end
    
end

