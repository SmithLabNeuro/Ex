function smoothedBins = exponentialSmoother(currTrialBins, alpha)
% Assumes currTrialBins is of num_dims x num_time_bins
    smoothedBins = nan(size(currTrialBins));
    for i =1:size(smoothedBins,2)
        % Set first value to currTrialBins first value
        if i == 1
            smoothedBins(:,i) = currTrialBins(:,i);
        else
            smoothedBins(:,i) = alpha*smoothedBins(:,i-1) + (1-alpha)*currTrialBins(:,i);
        end
    end
end