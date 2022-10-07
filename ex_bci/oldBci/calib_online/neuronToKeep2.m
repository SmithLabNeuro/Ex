function [neuronstokeep, modDepth] = neuronToKeep2(dat,varargin)

    %% make sure important stuff in path
    addpath(genpath('../structbuilders'))

    %% optional parameters
    nDelayBins = 0;
    
    goodRatesFlag = true;
    minRate = 1;
    
    goodDepthFlag = true;
    sigThresh = .8;
    nPermutations = 1000;
    binsize = 0.05;
    
    assignopts(who,varargin);

    %% get some important variables
    nNeurons = size(dat(1).counts,1);
    nTrials = length(dat);
    minBins = min([dat.nBins]);
    nTimePoints = minBins - 2*nDelayBins;
    targAngles = [dat.angle]';

    %% compute counts
    allCounts = nan(nNeurons,nTimePoints,nTrials);
    startIdx = nDelayBins+1;
    endIdx = nDelayBins+nTimePoints;
    for i_trial = 1:nTrials
        allCounts(:,:,i_trial) = dat(i_trial).counts(:,startIdx:endIdx);
    end
    allCounts = permute(allCounts,[1 3 2]);
    
    % get avg delay period FR for each neuron
    trialAvgCounts = squeeze(mean(allCounts,3));

    %% compute tuning curves
    [tuneCurves,~,targAngOrder] = getTuningCurves(trialAvgCounts,targAngles);

    %% compute modulation depth
    [~,~,~,modDepth,~] = fitCosineTuning(tuneCurves,targAngOrder);
    
    %% compute permutation test for modulation depth
    perm_modDepth = nan(nNeurons,nPermutations);
    for i_perm = 1:nPermutations
        if mod(i_perm,floor(nPermutations/10))==0
            fprintf('   Permutation test sample %d of %d\n',i_perm,nPermutations);
        end
        rngIdx = randperm(nTrials);
        [tmpCurves,~,tmpOrder] = getTuningCurves(trialAvgCounts,targAngles(rngIdx));
        [~,~,~,tmp,~] = fitCosineTuning(tmpCurves,tmpOrder);
        perm_modDepth(:,i_perm) = tmp;
    end
    
    %% decide on good neurons
    if goodRatesFlag
        goodRates = mean(trialAvgCounts,2)>(minRate*binsize);
    else
        goodRates = true(nNeurons,1);
    end
    
    if goodDepthFlag
        sigModDepth = sum(repmat(modDepth,1,nPermutations)>perm_modDepth,2)>floor(nPermutations*sigThresh);
    else
        sigModDepth = true(nNeurons,1);
    end
   
    neuronstokeep = goodRates & sigModDepth;
    
end