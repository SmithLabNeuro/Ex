function [ params, trainError ] = train_nb( trainX, trainY, varargin )
%
%   Input:
%       trainX   - features (nSamples x xDim)
%       trainY   - class labels (nSamples x 1)
%       distType - distribution of class-conditional densities
%                  'diagGauss' (default),'poisson', 'fullGauss',
%                  'sharedGauss'
%
%   Output:
%       params     - model parameters and distType
%                       - classPi, mu, and sigma for the gauss models
%                       - classPi, lambda for the poisson models
%       LL         - log likelhood
%       trainError - training error (as a percentage)
%
    
    distType = 'diaggauss';
    minVar = 0.01;
    priorType = 'set';
    assignopts(who,varargin);
    
    classLabels = unique(trainY);
    nClasses = length(classLabels);
    params.classLabels = classLabels;
    params.nClasses = nClasses;
    params.minVar = minVar;
    
    params.distType = distType;
    params.priorType = priorType;
    
    % keep track of sufficient statistics
    for ii=1:nClasses
        curr_class = trainX(trainY==classLabels(ii),:);
        params.suffStat(ii).sumX = sum(curr_class,1);
        params.suffStat(ii).N = size(curr_class,1);
        params.suffStat(ii).sumX2 = sum(curr_class.^2,1);
    end
    
    % fit class probabilities (classPi)
    switch lower(priorType)
        case 'set'
            for ii=1:nClasses
                params.classPi(ii) = sum(trainY==classLabels(ii))./length(trainY);
            end
        case 'equal'
            for ii=1:nClasses
                params.classPi(ii) = 1./nClasses;
            end
        otherwise
            error('Please specify a valid prior type')
    end        
       
    % fit model parameters    
    switch lower(distType)
        case 'diaggauss' % each class gets its own cov, but it's diagonal
            % fit means
            for ii=1:nClasses
                curr_class = trainX(trainY==classLabels(ii),:);
                params.Mu(ii,:) = mean(curr_class,1);
            end
            % fit covariance
            for ii=1:nClasses
                curr_class = trainX(trainY==classLabels(ii),:);
                params.Sigma(ii,:,:) = diag(max(diag(cov(curr_class,1)),minVar));
            end
        case 'shareddiag' % each class shares same cov, and it's diagonal
            % fit means
            for ii=1:nClasses
                curr_class = trainX(trainY==classLabels(ii),:);
                params.Mu(ii,:) = mean(curr_class,1);
            end
            % fit covariance
            for ii=1:nClasses
                curr_class = trainX(trainY==classLabels(ii),:);
                params.Sigma(ii,:,:) = diag(max(diag(cov(curr_class,1)),minVar));
            end
            newCov = diag(diag(squeeze(mean(params.Sigma,1))));
            for ii = 1:nClasses
                params.Sigma(ii,:,:) = newCov;
            end
        case 'fullgauss' % each class gets its own cov
            % fit means
            for ii=1:nClasses
                curr_class = trainX(trainY==classLabels(ii),:);
                params.Mu(ii,:) = mean(curr_class,1);
            end
            % fit covariance
            for ii=1:nClasses
                curr_class = trainX(trainY==classLabels(ii),:);
                params.Sigma(ii,:,:) = cov(curr_class,1);
            end
        case 'sharedgauss' % gaussian w/ same diagonal cov for each class
            % fit means
            for ii=1:nClasses
                curr_class = trainX(trainY==classLabels(ii),:);
                params.Mu(ii,:) = mean(curr_class,1);
            end
            % fit covariance
            modelCov = zeros(size(trainX,2));
            for ii=1:nClasses
                curr_class = trainX(trainY==classLabels(ii),:);
                modelCov = modelCov + params.classPi(ii) .* cov(curr_class,1);
            end
            for ii=1:nClasses
                % diagonalize the matrix and set a minimum variance
                params.Sigma(ii,:,:) = modelCov;
            end
        case 'poisson'
            % fit lambda
            for ii=1:nClasses
                curr_class = trainX(trainY==classLabels(ii),:);
                params.Lambda(ii,:) = mean(curr_class,1);
            end
            % set a minimum lambda
            params.Lambda = max(params.Lambda,minVar);
        otherwise
            error('Please specify a valid distribution type')
    end
    
    % calculate training error
    predictions = predict_nb(trainX,params);
    trainError = sum(predictions~=trainY)./length(trainY);
    
end

