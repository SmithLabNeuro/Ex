function [ params ] = train_hmm( trainX, trainY, varargin )
% trainX: xDim x nSamples x nTimePoints
% trainY: nSamples x 1
    

    stateLabels = unique(trainY); % i.e. classLabels
    params.stateLabels = stateLabels;
    nStates = length(params.stateLabels);
    params.nStates = nStates;
    
    stateTransMat = ones(nStates)./nStates;
    covType = 'full';
    assignopts(who,varargin);
    
    
    % probability for starting in a given state (assume flat prior)
    params.pi = ones(params.nStates,1)./params.nStates;
    
    % state transition matrix 
    params.T = stateTransMat;
    if sum(sum(params.T,2)~=1)>0
        error('Transition probabilities are not valid distributions (rows should sum to 1');
    end
    
    % fit gaussian model for each class
    for ii = 1:nStates
        currDat = trainX(:,trainY==stateLabels(ii),:);
        currDat = reshape(currDat,size(currDat,1),[]);
        params.Mu(ii,:) = mean(currDat,2);
        params.Sigma(ii,:,:) = cov(currDat',1);
    end
    
    % adjust to the specified covType
    switch lower(covType)
        case 'shareddiag'
            newCov = diag(diag(squeeze(mean(params.Sigma,1))));
            for ii = 1:nStates
                params.Sigma(ii,:,:) = newCov;
            end
        case 'diag'
            for ii = 1:nStates
                currCov = squeeze(params.Sigma(ii,:,:));
                params.Sigma(ii,:,:) = diag(diag(currCov));
            end
        case 'sharedfull'
            newCov = squeeze(mean(params.Sigma,1));
            for ii = 1:nStates
                params.Sigma(ii,:,:) = newCov;
            end
        case 'full'
        otherwise
            error('Please specify a valid covariance type (shareddiag, diag, sharedfull, full)');
    end
    
end
