function [ goodNeurons, modelPreds ] = train_lowdBayes( trainDat, goodNeurons, varargin )

    %% set optional arguments
    nDelayBins = 0;
    zDim = 2;
    covType = 'fullgauss';
    whichArray = 'both';
    nFolds = 15;
    useLDA = false;
    assignopts(who,varargin);
    
    
    %% get some important variables
    nNeurons = size(trainDat(1).counts,1);
    nTrials = length(trainDat);
    minBins = min([trainDat.nBins]);
    nTimePoints = minBins - nDelayBins;
    targAngles = [trainDat.angle]';
    
    
    %% determine which goodNeurons to keep based on whichArray specification
    switch lower(whichArray)
        case 'left'
            goodNeurons(1:96) = false;
        case 'right'
            goodNeurons(97:end) = false;
        case 'both'
        otherwise
            error('Specify a valid value for whichArray ("both","left","right")');
    end
    
    %% parse trainDat into a matrix of nNeurons x nTrials x nTimePoints
    allCounts = nan(nNeurons,nTimePoints,nTrials);
    startIdx = nDelayBins+1;
    endIdx = nDelayBins+nTimePoints;
    for i_trial = 1:nTrials
        allCounts(:,:,i_trial) = trainDat(i_trial).counts(:,startIdx:endIdx);
    end
    allCounts = permute(allCounts,[1 3 2]);
    
    % only keep good neurons
    allCounts = allCounts(goodNeurons,:,:);
    nNeurons = size(allCounts,1);
    
    % get avg counts on each trial
    trialCounts = squeeze(mean(allCounts,3));
    
    
    %% train LDA and project data into LDA space, compute training predictions
    if useLDA
        indCounts = reshape(allCounts,nNeurons,[])';
        indCountLabels = reshape(repmat(targAngles,1,size(allCounts,3)),1,[])';
        ldaParams = train_lda(indCounts,indCountLabels);
        ldaProjMat = ldaParams.projMat(:,1:zDim);
    else
        ldaProjMat = eye(nNeurons);
    end
    
    % project to LDA space
    lowd_x = ldaProjMat'*trialCounts;
    lowd_x = lowd_x';
    
    % compute training predictions
    train_params = train_nb(lowd_x,targAngles,'distType',covType,'priorType','equal');
    modelPreds.y_training = predict_nb(lowd_x,train_params)';
    modelPreds.y_true = targAngles';
    
    
    %% cross-validation
    fdiv = floor(linspace(1,nTrials+1,nFolds+1));
    yhat = nan(1,nTrials);
    
    for cvf = 1:nFolds
        testMask = false(1,nTrials);
        testMask(fdiv(cvf):(fdiv(cvf+1)-1)) = true;
        trainMask = ~testMask;
        trainX = trialCounts(:,trainMask);
        trainY = targAngles(trainMask);
        testX = trialCounts(:,testMask);
        
        % fit LDA
        if useLDA
            indCounts = reshape(trainX,nNeurons,[])';
            indCountLabels = reshape(repmat(trainY,1,size(trainX,3)),1,[])';
            ldaParams = train_lda(indCounts,indCountLabels);
            ldaProjMat = ldaParams.projMat(:,1:zDim);
        else
            ldaProjMat = eye(nNeurons);
        end
        lowd_trainX = ldaProjMat'*trainX; 
        lowd_trainX = lowd_trainX';
        lowd_testX = ldaProjMat'*testX; 
        lowd_testX = lowd_testX';
        
        cvf_params = train_nb(lowd_trainX,trainY,'distType',covType,'priorType','equal');
        yhat(testMask) = predict_nb(lowd_testX,cvf_params);
    end
    modelPreds.y_cv = yhat;
    
%     %% train LDA and project data into LDA space
%     indCounts = reshape(allCounts,nNeurons,[])';
%     indCountLabels = reshape(repmat(targAngles,1,size(allCounts,3)),1,[])';
%     ldaParams = train_lda(indCounts,indCountLabels);
%     ldaProjMat = ldaParams.projMat(:,1:zDim);
%     
%     % sum over timepoints
%     allCounts = squeeze(sum(allCounts,3));
%     
%     % project to LDA space
%     lowd_x = ldaProjMat'*allCounts;
%     
%     %% train bayes classifier on LDA projections
%     [nb_params,~,heldOutPred] = crossvalidate_nb(lowd_x',targAngles,...
%         'nFolds',15,'distType',covType,'priorType','equal','useRandomSeed',true);
    
    
    %% evaluate cross-validated prediction accuracy
    predAcc = sum(modelPreds.y_true==modelPreds.y_cv)./length(modelPreds.y_true).*100;
    fprintf('Cross-validated prediction accuracy for %s array(s): %.2f\n',whichArray,predAcc);
%     
%     modelPreds.y_true = heldOutPred.y;
%     modelPreds.y_cv = heldOutPred.yhat;
%     modelPreds.y_training = predict_nb(heldOutPred.X,nb_params);
    
end

