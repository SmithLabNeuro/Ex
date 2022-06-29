function [ predictions, diff_logP, log_post ] = predict_nb( testX, params )
%
%   Input:
%       testX - features (nSamples x xDim)
%       params - trained parameters
%
%   Output:
%       predictions - (nSamples x 1)
%

    nClasses = params.nClasses;
    distType = params.distType;
    classLabels = params.classLabels;
    
    [nSamples,nDims] = size(testX);
    logP = zeros(nSamples,nClasses);
    
    % un-normalized posterior probability
    if strcmpi(distType,'poisson')
        for ii=1:nClasses
            Lambda = params.Lambda(ii,:);
            classPi = params.classPi(ii);
            logP(:,ii) = testX*log(Lambda') - sum(Lambda) + log(classPi);
        end
    else % gaussian
        for ii=1:nClasses
            Mu = params.Mu(ii,:);
            Sigma = squeeze(params.Sigma(ii,:,:));
            classPi = params.classPi(ii);
            centeredX = bsxfun(@minus,testX,Mu);
            part1 = -1/2 * diag(centeredX/Sigma*centeredX');
            part2 = -1/2 * log(det(Sigma));
            part3 = -nDims/2 * log(2*pi);
            part4 = log(classPi);
            logP(:,ii) = part1 + part2 + part3 + part4;
        end
    end
    
    % log likelihood
    consts = max(logP, [], 2);
    normLogL = bsxfun(@minus, logP, consts);
    p1_plus_p2 = log(sum(exp(normLogL), 2)) + consts;
    
    % log posterior probability
    log_post = bsxfun(@minus,logP,p1_plus_p2);
    
    [~, maxIdx] = max(logP,[],2);
    predictions = classLabels(maxIdx);
    diff_logP = logP(:,2)-logP(:,1);
    
    
    
end

