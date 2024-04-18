function [params] = fit_LDA(trainX, trainY, numDim)
% Fits a linear discriminant model to provided trainX and trainY
% find projection vector and project X 
    % This backslash is used in place of inv to account for singular Sw
    % matrices. If singular, will use the least squares estimate intead of
    % throwing an error which inv() will do.
    discMdl = fitcdiscr(trainX, trainY); 
    sB = discMdl.BetweenSigma; 
    sW = discMdl.Sigma; 
    [V, eigVls] = eig(sB, sW); 
    % these are in reverse order for some reason... let's fix that 
    [eigVls, order] = sort(diag(eigVls),'descend'); 
    V = V(:,order); % find the rank (# labels - 1 or #dims) and extract w accordingly 
    params.projMat = V;
    params.projVec = V(:,1:numDim);
    params.projData = trainX*params.projVec;
    params.sB = sB;
    params.sW = sW;
end