function smoothedBins = exponentialSmoother(currTrialBins, alpha, initialValue)
% Assumes currTrialBins is of num_dims x num_time_bins
% InitialValue corresponds to a "seed" value can that bias future smoothed
% values.
    smoothedBins = nan(size(currTrialBins));
    for i =1:size(smoothedBins,2)
        % Set first value to provided initial value instead of currTrialBins first value
        if i == 1
            smoothedBins(:,i) = initialValue;
        else
            smoothedBins(:,i) = alpha*smoothedBins(:,i-1) + (1-alpha)*currTrialBins(:,i);
        end
    end
end