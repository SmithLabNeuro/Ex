function [trimmedDat, channelsKeep] = preprocessDatWithNoStimTrial(datStruct, nevLabelledData, channelNumbersUse, binSizeMs, frThresh, ffThresh, coincThresh, coincidenceTimeMs)
% this function is responsible for:
% - removing low firing rate channels
% - removing high fano factor channels
% - computing coincident channels using *no uStim* calibration data and removing
%   them

% extract no uStim trials
% TODO: write using the regular expression
stimChan = cellfun(@(dS) regexp(dS, 'xippmexStimChan*=((\d*\s*)*);.*', 'tokens'), {datStruct.text});
stimChan = cellfun(@(a) a{1}, stimChan, 'UniformOutput', false);
datStruct = datStruct(strcmp(stimChan, '0')|strcmp(stimChan, '0  0')|strcmp(stimChan, '0  0  0')|strcmp(stimChan, '0  0  0  0')|strcmp(stimChan, '0  0  0  0  0'));

samplingRate = 30e3; % Ripple Hz
spikeTimes = cellfun(@(firstSpike, spikeTimeDiffs) [firstSpike; firstSpike+cumsum(uint64(spikeTimeDiffs))], {datStruct.firstspike}, {datStruct.spiketimesdiff}, 'uni', 0);
spikeChannels = cellfun(@(spikeInfo) spikeInfo(:, 1), {datStruct.spikeinfo}, 'uni', 0);
trialStartEndsSec = cat(1,datStruct.time);
trialLengthsSec = diff(trialStartEndsSec, [], 2);
totalTrialTime = sum(trialLengthsSec);

channelNumbers = unique(cat(1, spikeChannels{:}));
% grab only the channels from the array that actually had a signal
channelNumbersUse = channelNumbersUse(ismember(channelNumbersUse, channelNumbers));
channelNumbersNotUse = channelNumbers(~ismember(channelNumbers, channelNumbersUse))';
numChannels = length(channelNumbersUse);
msToS = 1000;
coincidenceTimeSamples = coincidenceTimeMs/msToS*samplingRate;

binSizeS = binSizeMs/msToS; % s; having the same bin size is important for Fano factor, not for firing rate
spikesPerChanPerTrial = cellfun(...
                @(trialSpikeTimes, chanNums)...
                histcounts2(trialSpikeTimes, chanNums, trialSpikeTimes(1):binSizeS*samplingRate:max(trialSpikeTimes)+1, 1:max(channelNumbers)+1),...
                spikeTimes, spikeChannels, 'uni', 0);
spikesPerChanPerTrialMat = cat(1, spikesPerChanPerTrial{:});
spikesPerChanPerTrialMat = spikesPerChanPerTrialMat(:, channelNumbersUse);
totalSpikesPerChan = sum(spikesPerChanPerTrialMat, 1);


chanFR = totalSpikesPerChan/totalTrialTime;
chanFF = var(spikesPerChanPerTrialMat, [], 1)./mean(spikesPerChanPerTrialMat, 1);

% finding coincidence--this is gonna be a channel-by-channel approach
coincCh = zeros(numChannels);
chNumSpikes = zeros(numChannels, 1);
spikeTimesMat = cat(1, spikeTimes{:});
spikeChannelsMat = cat(1, spikeChannels{:});
fprintf('\n')
fprintf('Processing %d channels for coincidence\n', numChannels);
for chInd = 1:numChannels
    ch = channelNumbersUse(chInd);
    if ~mod(chInd, 10)
        fprintf('  %d channels left\n', numChannels-chInd);
    end
    % chNumSpikes(chInd) = sum(nevLabelledData(:, 1)==ch); % spike + noise (NASNet output base)
    chNumSpikes(chInd) = sum(spikeChannelsMat==ch); % spike only (dat structure base)

    for ch2Ind = chInd+1:numChannels
        % this is an order of magnitude faster than trying to find all the
        % pairwise differences between all spikes (perhaps unsurprisingly),
        % and semi-relies on the coincidence times being < a waveform
        % length, so this will never be true between two waveforms on one
        % channel
        ch2 = channelNumbersUse(ch2Ind);
        % dfsBothCh = diff(samplingRate*nevLabelledData(nevLabelledData(:,1)==ch | nevLabelledData(:, 1)==ch2, 3))<=coincidenceTimeSamples;
        dfsBothCh = diff(spikeTimesMat(spikeChannelsMat==ch | spikeChannelsMat==ch2))<=coincidenceTimeSamples;
        numCloseSpks = sum(dfsBothCh);
        
        coincCh(chInd, ch2Ind) = coincCh(chInd, ch2Ind) + numCloseSpks;
        coincCh(ch2Ind, chInd) = coincCh(chInd, ch2Ind);
    end
end
fprintf('Coincidence measurement done\n');

coincProp = coincCh./chNumSpikes;
% below, each entry is the number of channels for which the current
% channels has more that coincThresh fraction spikes coincident
numChCoincidentTo = sum(coincProp>coincThresh, 2); 
% below, each entry is the number of channels that had more than
% coincThresh fraction spikes coincident with this channel's spikes
numChCoincidentWith = sum(coincProp>coincThresh, 1);
% so that we can keep the original coincProp around (for potential later
% use but also debugging purposes...)
coincPropNew = coincProp;

chansRemaining = channelNumbersUse;
while sum(numChCoincidentTo)>0 && sum(numChCoincidentWith)>0
    % compute the channel that is either maximally coincident to other
    % channels or with which a maximal number of channels are coincident
    [mxCoinc, mxCoincidenceCh] = max([numChCoincidentTo; numChCoincidentWith']);
    % what happens if there are multiple channels with the max? We try to
    % find the one whose total is maxed out
    if (sum(numChCoincidentTo == mxCoinc) + sum(numChCoincidentWith == mxCoinc)) > 1
        sumCoinc = numChCoincidentTo + numChCoincidentWith';
        [~, mxCoincidenceCh] = max(sumCoinc);
        % now we bias ourselves to just grabbing the first of them--I'd do
        % random but I'd need to be careful about having a set seed every
        % time, so I think this is fine
    end
    % in case it's the numChCoincidentWith, gotta mod it back to the range
    % of channels
    mxCoincidenceCh = mod(mxCoincidenceCh, length(numChCoincidentTo));
    % Matlab's 1-indexed, mod returns zero for the length matching, so fix
    % it here
    if mxCoincidenceCh == 0
        mxCoincidenceCh = length(numChCoincidentTo);
    end
    
    % remove said channel from remaining channels
	chansRemaining(mxCoincidenceCh) = [];
    
    % gotta recompute numChCoincidentTo and numChCoincidentWith because
    % different entries might have included this channel or not, so we go
    % back to coincProp and remove the channel there before recomputing
    coincPropNew(:, mxCoincidenceCh) = [];
    coincPropNew(mxCoincidenceCh, :) = [];
    
    % below, each entry is the number of channels for which the current
    % channels has more that coincThresh fraction spikes coincident
    numChCoincidentTo = sum(coincPropNew>coincThresh, 2);
    % below, each entry is the number of channels that had more than
    % coincThresh fraction spikes coincident with this channel's spikes
    numChCoincidentWith = sum(coincPropNew>coincThresh, 1); 
    
end

% now we can remove any remaining firing rate and fano factor issue
% channels
chanRemoveLogical = chanFR < frThresh | chanFF > ffThresh | ~ismember(channelNumbersUse, chansRemaining);
channelsRemove = channelNumbersUse(chanRemoveLogical);
channelsKeep = channelNumbersUse(~chanRemoveLogical);

channelsRemove = [channelsRemove channelNumbersNotUse];

nevLabelledData(ismember(nevLabelledData(:, 1), channelsRemove), :) = [];
trimmedDat = nev2dat(nevLabelledData, 'nevreadflag', true);


