function [countsCurrentBin, countsNextBin] = countBinnedSpikesPerChannel(timestampsThresholdCrossings, channelsToCount, timePtBinStart, samplesPerBin, nasNetParams, waveformsThresholdCrossings, gamma)
% computes number of spikes in current and next bin
% input:
%   timestampsThresholdCrossings : channel x 1 cell array of timestamps for
%       threshold crossings, where each cell contains all timestamps for
%       that channel
%   channelsToCount : channels x 1 logical array where a true value
%       indicates a channel whose bins to count
%   timePtBinStart : scalar indicating timestamp for the start of the bin
%   samplesPerBin : scalar indicating number of samples in a bin
%   nasNetwork : a filepath indicating the NAS network to use for filtering
%      out waveforms; if empty or this and following arguments aren't
%      input, then NAS network is not used for filtering and the following
%      inputs are ignored
%   waveformsThresholdCrossings : channel x 1 cell array of waveforms for
%       threshold crossings, where each cell contains a (num timestamps x
%       samples per waveform) array of waveforms associated with each
%       timestamp for that channel in timestampsThresholdCrossings
%   gamma : gamma value to use for NAS network to allow waveforms through
% output: 
%   countsCurrentBin : channel x 1 array of spike counts in the current bin
%   countsNextBin : channel x 1 array of spike counts in the next bin
% global tmstpNasSpkAll
if nargin<5 || isempty(nasNetParams)
    useNas = false;
else
    useNas = true;
    [w1, b1, w2, b2] = deal(nasNetParams.w1, nasNetParams.b1, nasNetParams.w2, nasNetParams.b2);
end

% only use data from channels that were good for calibration
tmstpGoodChSpk = timestampsThresholdCrossings(channelsToCount);

if useNas
    % run waveforms through NAS net
    waveformsThresholdCrossings = waveformsThresholdCrossings(channelsToCount);
    tmstpFiltered = cellfun(@(wvForms, tms) tms(runNASNetContinuous(w1, b1, w2, b2, wvForms, gamma)), waveformsThresholdCrossings, tmstpGoodChSpk, 'uni', 0);
else
    % pass all waveeforms through
    tmstpFiltered = tmstpGoodChSpk;
end
% tmstpNasSpkAll = [tmstpNasSpkAll tmstpFiltered];

% grab timestamps in current bin
spikesThisBinByChannel = cellfun(@(x) x>=timePtBinStart & x<timePtBinStart+samplesPerBin, tmstpFiltered, 'uni', 0);
% grab timestamps in next bin
spikesNextBinByChannel = cellfun(@(x) x>=timePtBinStart+samplesPerBin, tmstpFiltered, 'uni', 0);

% compute spike counts for current and next bin
countsCurrentBin = countMatFromTimestamps(spikesThisBinByChannel);
countsNextBin = countMatFromTimestamps(spikesNextBinByChannel);


function channelCountMat = countMatFromTimestamps(cellOfTimestamps)
% subfunction to give a channel x 1 array of counts given a channel x 1
% cell array of timestamps (where each channel's cell contains all
% timestamps for that channel)
countsPerChannelCell = cellfun(@(x) sum(x, 2), cellOfTimestamps, 'uni', 0);
countsExistChannel = ~cellfun('isempty', countsPerChannelCell);
channelCountMat = zeros(length(countsPerChannelCell), 1);
channelCountMat(countsExistChannel) = [countsPerChannelCell{countsExistChannel}];