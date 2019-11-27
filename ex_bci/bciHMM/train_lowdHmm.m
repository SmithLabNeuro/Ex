function [ modelParams ] = train_lowdHmm( trainDat, modelParams, varargin )
    

    %% set optional arguments
    nDelayBins = 0;
    zDim = 2;
    hmmCovType = 'sharedfull';
    evalTrainingErr = false;
    visualizeLDA = false;
    useLDA = false;
    assignopts(who,varargin);
    
    modelParams.nDelayBins = nDelayBins;
    modelParams.zDim = zDim;
    
    
    %% get some important variables
    nNeurons = size(trainDat(1).counts,1);
    nTrials = length(trainDat);
    minBins = min([trainDat.nBins]);
    nTimePoints = minBins - 2*nDelayBins;
    targAngles = [trainDat.angle]';
    unTargAngles = unique(targAngles);
    
    
    %% parse trainDat into a matrix of nNeurons x nTrials x nTimePoints
    allCounts = nan(nNeurons,nTimePoints,nTrials);
    startIdx = nDelayBins+1;
    endIdx = nDelayBins+nTimePoints;
    for i_trial = 1:nTrials
        allCounts(:,:,i_trial) = trainDat(i_trial).counts(:,startIdx:endIdx);
    end
    allCounts = permute(allCounts,[1 3 2]);
    
    % only keep good neurons
    allCounts = allCounts(modelParams.goodneuron,:,:);
    nNeurons = size(allCounts,1);
    
    
    %% compute centering mean and angle means; center allCounts
    trialMeans = squeeze(nanmean(allCounts,3)); % nNeurons x nTrials
    calib_angMeans = nan(nNeurons,length(unTargAngles));
    for i_ang = 1:length(unTargAngles)
        currDat = trialMeans(:,targAngles==unTargAngles(i_ang));
        calib_angMeans(:,i_ang) = mean(currDat,2);
    end
    centeringMean = mean(calib_angMeans,2);
    allCounts = bsxfun(@minus,allCounts,centeringMean);
    modelParams.calib_angMeans = calib_angMeans;
    modelParams.centeringMean = centeringMean;
    modelParams.targetAngles = unTargAngles;
    
    
    %% train LDA and project data into LDA space
    if useLDA
        indCounts = reshape(allCounts,nNeurons,[])';
        indCountLabels = reshape(repmat(targAngles,1,size(allCounts,3)),1,[])';
        ldaParams = train_lda(indCounts,indCountLabels);
        %ldaParams = train_lda(squeeze(mean(allCounts,3))',targAngles);
        modelParams.ldaProjMat = ldaParams.projMat(:,1:zDim);
    else
        modelParams.ldaProjMat = eye(nNeurons);
    end
    
    lowd_x = nan(size(modelParams.ldaProjMat,2),nTrials,nTimePoints);
    for i_t = 1:nTimePoints
        lowd_x(:,:,i_t) = modelParams.ldaProjMat'*squeeze(allCounts(:,:,i_t));
    end
    
    
    %% train hmm
    transMatrix = eye(length(unique(targAngles)));
    modelParams.hmmParams = train_hmm(lowd_x,targAngles,'stateTransMat',transMatrix,'covType',hmmCovType);
    
    
    %% evaluate training error of model over time
    if evalTrainingErr
        for i_trial = 1:nTrials
            [predAngle(i_trial,:),~] = hmm_filter(squeeze(lowd_x(:,i_trial,:)),modelParams.hmmParams);
        end
        for i_t = 1:nTimePoints
            trainAcc(i_t) = sum(predAngle(:,i_t)==targAngles)./length(targAngles);
        end
        figure; 
        plot(startIdx:endIdx,trainAcc,'o-');
        xlabel('bin from target offset'); ylabel('Training accuracy'); 
        title('LDA+HMM');
    end
    
    
    %% makes plots to visualize training data
    if visualizeLDA
        figure; pos=get(gcf,'Position'); set(gcf,'Position',pos.*[1 1 2 2]);
        
        % plot hmm means
        subplot(2,2,1); hold on;
        for i_ang = 1:modelParams.hmmParams.nStates
            meanVal = modelParams.hmmParams.Mu(i_ang,:);
            plot(meanVal(1),meanVal(2),'.','MarkerSize',15);
        end
        legend(cellstr(num2str(modelParams.hmmParams.stateLabels)),'Location','Best');
        title('HMM means'); xlabel('LDA_1'); ylabel('LDA_2');
        
        % plot each trial
        subplot(2,2,2); hold on;
        lowd_x_avg = squeeze(mean(lowd_x,3));
        for i_ang = 1:modelParams.hmmParams.nStates
            currDat = lowd_x_avg(:,modelParams.hmmParams.stateLabels(i_ang)==targAngles);
            plot(currDat(1,:),currDat(2,:),'.');
        end
        title('All calibration trials'); xlabel('LDA_1'); ylabel('LDA_2');
        
        % plot each of the 4 targets for BCI
        subplot(2,2,3); hold on;
        bciAngs = modelParams.hmmParams.stateLabels(1:4:end);
        for i_ang = 1:length(bciAngs)
            currDat = lowd_x_avg(:,bciAngs(i_ang)==targAngles);
            plot(currDat(1,:),currDat(2,:),'.');
        end
        legend(cellstr(num2str(bciAngs)),'Location','Best');
        title('Calibration trials for some targets'); xlabel('LDA_1'); ylabel('LDA_2');
    end
    
    
end

