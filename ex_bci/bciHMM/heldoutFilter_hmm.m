function [ progMet, postProbs ] = heldoutFilter_hmm( X, Y, varargin )
%
% trainX: xDim x nSamples x nTimePoints
% trainY: nSamples x 1
% postProbs: nStates x nTimePoints x nSamples
% progMet: nSamples x nTimePoints
%


    % some preliminary and optional parameters
    stateLabels = unique(Y); % i.e. classLabels
    nStates = length(stateLabels);
    [~,nSamples,nT] = size(X);
    nFolds = nSamples/nStates;
    if mod(nFolds,1)~=0
        error('Number of folds is not an integer')
    end
    stateTransMat = ones(nStates)./nStates;
    covType = 'full';
    assignopts(who,varargin);
    
    % determine cross-validation folds
    fdiv = floor(linspace(1,nSamples+1,nFolds+1));
    
    postProbs = nan(nStates,nT,nSamples);
    progMet = nan(nSamples,nT);
    for cvf=1:nFolds
        % separate train and test data
        testMask = false(1,nSamples);
        testMask(fdiv(cvf):(fdiv(cvf+1)-1)) = true;
        trainMask = ~testMask;
        trainX = X(:,trainMask,:);
        trainY = Y(trainMask);
        testX = X(:,testMask,:);
        testY = Y(testMask);
        
        % fit model to training data
        cvf_params = train_hmm(trainX,trainY,'stateTransMat',stateTransMat,'covType',covType);
        
        % filter the test data using trained model
        tmpPost = nan(nStates,nT,size(testX,2));
        tmpProg = nan(size(testX,2),nT);
        for ii = 1:length(testY)
              [~,tmpPost(:,:,ii)] = hmm_filter(squeeze(testX(:,ii,:)),cvf_params);
              tmpProg(ii,:) = bciProgMetric(tmpPost(:,:,ii)',stateLabels,testY(ii));
        end
        progMet(testMask,:) = tmpProg;
        postProbs(:,:,testMask) = tmpPost;
    end
    
    
    
end

