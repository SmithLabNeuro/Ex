function [binSpikeCountOverall, binSpikeCountNextOverall] = countBinnedSpikesOverall(binSpikeCountOverall, binSpikeCountNextOverall, tmstpInit, goodChannelInds, timePtBinStart, samplesPerBin, nasNetParams, prevWaveforms, gamma)

[countsPerChannel, countsPerChannelNextBin] = countBinnedSpikesPerChannel(tmstpInit, goodChannelInds, timePtBinStart, samplesPerBin, nasNetParams, prevWaveforms, gamma);
binSpikeCountOverall = binSpikeCountOverall + countsPerChannel;
binSpikeCountNextOverall = binSpikeCountNextOverall + countsPerChannelNextBin;
