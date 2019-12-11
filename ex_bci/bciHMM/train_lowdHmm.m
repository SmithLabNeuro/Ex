function [ modelParams ] = train_lowdHmm( trainDat, modelParams, varargin )
    

    %% set optional arguments
    nDelayBins = 8;
    zDim = 5;
    hmmCovType = 'shareddiag';
    evalTrainingErr = false;
    transMatType = 'diag';
    diag_weight = 0.9;
    visualizeLDA = false;
    assignopts(who,varargin);
    
    modelParams.nDelayBins = nDelayBins;
    modelParams.zDim = zDim;
    
    
    %% get some important variables
    nNeurons = size(trainDat(1).counts,1);
    nTrials = length(trainDat);
    minBins = min([trainDat.nBins]);
    nTimePoints = minBins - nDelayBins;
    targAngles = [trainDat.angle]';
    modelParams.allangles = unique(targAngles);
    
    
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
    modelParams.centeringMean = squeeze(mean(mean(allCounts,2),3));
    
    %% train LDA and project data into LDA space
    indCounts = reshape(allCounts,nNeurons,[])';
    indCountLabels = reshape(repmat(targAngles,1,size(allCounts,3)),1,[])';
    ldaParams = train_lda(indCounts,indCountLabels);
    %ldaParams = train_lda(squeeze(mean(allCounts,3))',targAngles);
    modelParams.ldaProjMat = ldaParams.projMat(:,1:zDim);
    
    lowd_x = nan(zDim,nTrials,nTimePoints);
    for i_t = 1:nTimePoints
        lowd_x(:,:,i_t) = modelParams.ldaProjMat'*squeeze(allCounts(:,:,i_t));
    end
    
    
    %% train hmm
    if strcmpi(transMatType,'diag')
        transitionMatrix = eye(length(unique(targAngles)));
    elseif strcmpi(transMatType,'offDiag')
        transitionMatrix = zeros(length(unique(targAngles)));
        offdiag_weight = (1-diag_weight)/2;
        for ii = 1:length(unique(targAngles))
            transitionMatrix(ii,ii) = diag_weight;
            if ii==1
                transitionMatrix(ii,end) = offdiag_weight;
            else
                transitionMatrix(ii,ii-1) = offdiag_weight;
            end
            if ii==length(unique(targAngles))
                transitionMatrix(ii,1) = offdiag_weight;
            else
                transitionMatrix(ii,ii+1) = offdiag_weight;
            end
        end
    else
        error('Incorrect "transMatType" argument');
    end
    modelParams.hmmParams = train_hmm(lowd_x,targAngles,'stateTransMat',transitionMatrix,'covType',hmmCovType);
    
    
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

