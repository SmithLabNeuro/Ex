function [ postProb ] = hmm_onlineFilter( testX, prevPostProb, params )
%  testX        : nNeurons x 1 (a single timpoint within a trial)
%  prevPostProb : nStates x 1
%  params       : parameters fit using train_hmm
%  postProb     : nStates x 1

    %% important variables
    nStates = length(params.hmmParams.stateLabels);
    priorProb = params.hmmParams.pi;
    
    %% preprocess data (keep good neurons + project into LDA space)
    testX = testX - params.centeringMean;
    testX = params.ldaProjMat'*testX;
    
    %% forward algorithm
    alphaVals = nan(nStates,1);
    % if this is the first observation
    if isempty(prevPostProb)
        for s = 1:nStates
            emissionProb = mvnpdf(testX,params.hmmParams.Mu(s,:)',squeeze(params.hmmParams.Sigma(s,:,:)));
            alphaVals(s) = priorProb(s)*emissionProb;
        end
    else
        for s = 1:nStates
            emissionProb = mvnpdf(testX,params.hmmParams.Mu(s,:)',squeeze(params.hmmParams.Sigma(s,:,:)));
            transProbs = squeeze(params.hmmParams.T(:,s));
            alphaVals(s) = emissionProb*(transProbs'*prevPostProb);
        end
    end
    
    % scale alphas to get posterior probabilities
    scaleFact = sum(alphaVals);
    postProb = alphaVals./scaleFact;
    
end

