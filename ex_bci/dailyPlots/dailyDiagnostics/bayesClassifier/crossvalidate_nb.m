function [ params, cv_predError, heldOutPred ] = crossvalidate_nb( X, Y, varargin )
%CV_LDA2 Summary of this function goes here
%   Inputs:
%       X             - features (nSamples x xDim x nTimepoints)
%       Y             - class labels (nSamples x 1)
%       distType      - 'diagGauss' (default),'poisson', 'fullGauss', 'sharedGauss'
%       nFolds        - number of cross validation folds (optional, default is 10)
%       useRandomSeed - use a random number seed for randomizing samples (logical)
%       minVar        - minimum variance threshold
%       priorType     - use data proportions ('setprior'), or use flat priors ('equalprior')
%
%   Output:
%       params       - model parameters and distType
%                         - classPi, mu, and sigma for the gaussian models
%                         - classPi, lambda for the poisson models
%       cv_predError - prediction error across the n folds
%

    nFolds = 10;
    distType = 'diaggauss';
    useRandomSeed = false;
    minVar = 0.01;
    priorType = 'set';
    assignopts(who,varargin);
    
    nSamples = length(Y);
    
    % randomly reorder data points
    if useRandomSeed
        rng(0);
    end
    rng_idx = randperm(nSamples);
    X = X(rng_idx,:);
    Y = Y(rng_idx);
    fdiv = floor(linspace(1,nSamples+1,nFolds+1));
    
    % train on all data to get params
    [params,~] = train_nb(X,Y,'distType',distType,'minVar',minVar,'priorType',priorType);
    
    cv_predError = zeros(1,nFolds);
    heldOut = nan(2,nSamples);
    % cross-validation on nFolds
    for cvf=1:nFolds
        testMask = false(1,nSamples);
        testMask(fdiv(cvf):(fdiv(cvf+1)-1)) = true;
        trainMask = ~testMask;
        trainX = X(trainMask,:);
        trainY = Y(trainMask);
        testX = X(testMask,:);
        testY = Y(testMask);
        heldOut(1,fdiv(cvf):(fdiv(cvf+1)-1)) = testY;
        
        % fit model training data
        [cvf_params,~] = train_nb(trainX,trainY,'distType',distType,'minVar',minVar,'priorType',priorType);
        
        % evaluate model on test data
        test_predY = predict_nb(testX,cvf_params);
        heldOut(2,fdiv(cvf):(fdiv(cvf+1)-1)) = test_predY;
        cvf_PE = sum(test_predY~=testY)./length(testY);
        cv_predError(cvf) = cvf_PE;
    end
    
    heldOutPred.y = Y';
    heldOutPred.yhat = heldOut(2,:);
    heldOutPred.X = X;
    
end

