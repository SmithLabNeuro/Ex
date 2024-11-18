
function smoothedBins = exponentialSmoother(currTrialBins, alpha, initialValue)
% Assumes currTrialBins is of num_dims x num_time_bins
% InitialValue corresponds to a "seed" value can that bias future smoothed
% values.
    smoothedBins = nan(size(currTrialBins));
    for i =1:size(smoothedBins,2)
        % Set first value to provided initial value instead of currTrialBins first value
        if i == 1
            if isnan(initialValue) %no initial seed is provided use the current value
                smoothedBins(:,i) = currTrialBins(:,i);
            else
                smoothedBins(:,i) = initialValue;
            end
        else
            smoothedBins(:,i) = (1-alpha)*smoothedBins(:,i-1) + alpha*currTrialBins(:,i);
        end
    end
end