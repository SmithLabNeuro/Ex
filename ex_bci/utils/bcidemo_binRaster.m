function binned = bcidemo_binRaster(raster,bin_size,time_before,time_after)


% This function takes a raster (output of bcidemo_calculateRaster) and bins 
% the spikes into specified time intervals (bins) of a given size.
% 
% Parameters:
% 
    % raster: A cell array where each element is a vector containing spike times 
    % for a trial or unit. Each vector represents the times 
    % (in seconds) at which spikes occurred during the trial or 
    % for the specific unit.

    % bin_size: A scalar specifying the size of the bins in which the 
    % spikes will be grouped (seconds).
    
    % time_before: A scalar specifying the time before to use in the
    % binning (seconds).
  
    % time_after: A scalar specifying the time before to use in the binning
    %(seconds).

    % Returns:
    % binned: A matrix of size (length(edges)-1, length(raster)), where 
    % each column corresponds to a trial. Each row corresponds to the count 
    % of spikes within the respective time bin.

edges = time_before:bin_size:time_after;
assert((edges(end)-time_after)< 10^-10)
binned = nan(length(edges)-1,length(raster));

for i = 1:length(raster)

    % Calculate the counts within each bin
    counts = histcounts(raster{i}, edges);

    binned(:,i) = counts./bin_size;
end


end