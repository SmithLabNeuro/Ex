function [ x, y, RTs ] = binData( spikeStruct, binSize )
%BINDATA Summary of this function goes here
%   Detailed explanation goes here
    
    nNeurons = size(spikeStruct(1).counts_1ms,1);
    nTrials = length(spikeStruct);

    fixOffCode = 3;
    saccadeCode = 141;
    targOffCode = 100;
    
    RTs = nan(nTrials,1);
    delayPeriod = nan(nTrials,1);
    for i_trial = 1:nTrials
        spikeStruct(i_trial).nT = size(spikeStruct(i_trial).counts_1ms,2);
        spikeStruct(i_trial).targAngle = spikeStruct(i_trial).params.trial.angle;
        trialCodes = spikeStruct(i_trial).trialcodes;
        saccadeTime = trialCodes(ismember(trialCodes(:,2),saccadeCode),3);
        fixOffTime = trialCodes(ismember(trialCodes(:,2),fixOffCode),3);
        targOffTime = trialCodes(ismember(trialCodes(:,2),targOffCode),3);
        delayPeriod(i_trial) = fixOffTime(1) - targOffTime(1);
        if ~isempty(saccadeTime)
            RTs(i_trial) = saccadeTime - fixOffTime(1);
        end
    end
    
    maxT = max([spikeStruct.nT]);
    tmpSpikes = nan(nNeurons,maxT,nTrials);
    
    for i_trial=1:nTrials
        currT = spikeStruct(i_trial).nT;
        nBins(i_trial) = floor(currT/binSize);
        last_idx = nBins(i_trial)*binSize;
        tmpSpikes(:,1:last_idx,i_trial) = spikeStruct(i_trial).counts_1ms(:,1:last_idx);
    end

    maxBins = max(nBins);
    x = nan(nNeurons,maxBins,nTrials);
    for i_bin = 1:maxBins
        currIdx = (binSize*(i_bin-1)+1):(binSize*i_bin);
        x(:,i_bin,:) = sum(tmpSpikes(:,currIdx,:),2);
    end
    
    y = [spikeStruct.targAngle];
    
end

