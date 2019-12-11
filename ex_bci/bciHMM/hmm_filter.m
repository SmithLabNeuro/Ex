function [ pred, postProbs ] = hmm_filter( testX, params )
%  testX     : xDim x nTimePoints (a single sample)
%  params    : parameters fit using train_hmm
%  pred      : 1 x nTimePoints
%  postProbs : nStates x nTimePoints

    stateLabels = params.stateLabels;
    nStates = length(stateLabels);
    priorProb = params.pi;
    
    nT = size(testX,2);
    alphaVals = nan(nStates,nT);
    scaleFact = nan(1,nT);
    
    isDiagCov = isdiag(squeeze(params.Sigma(1,:,:)));
    
    % set initial alphaVals
    for s = 1:nStates
        if isDiagCov
            emissionProb = mvnpdf(testX(:,1),params.Mu(s,:)',diag(squeeze(params.Sigma(s,:,:)))');
        else
            emissionProb = mvnpdf(testX(:,1),params.Mu(s,:)',squeeze(params.Sigma(s,:,:)));
        end
        alphaVals(s,1) = priorProb(s)*emissionProb;
    end
    scaleFact(1) = sum(alphaVals(:,1));
    alphaVals(:,1) = alphaVals(:,1)./scaleFact(1);
    
    for t = 2:nT
        prevAlphas = alphaVals(:,t-1);
        currData = testX(:,t);
        for s = 1:nStates
            if isDiagCov
                emissionProb = mvnpdf(currData,params.Mu(s,:)',diag(squeeze(params.Sigma(s,:,:)))');
            else
                emissionProb = mvnpdf(currData,params.Mu(s,:)',squeeze(params.Sigma(s,:,:)));
            end
            transProbs = squeeze(params.T(:,s));
            alphaVals(s,t) = emissionProb*(transProbs'*prevAlphas);
        end
        scaleFact(t) = sum(alphaVals(:,t));
        alphaVals(:,t) = alphaVals(:,t)./scaleFact(t);
    end
    
    postProbs = alphaVals;
    [~,maxIdx] = max(postProbs,[],1);
    pred = stateLabels(maxIdx);
    
end

