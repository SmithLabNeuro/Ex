function [ params ] = train_lda( trainX, trainY )
    % trainX: nSamples x xDim
    % trainY: nSamples x 1

    % useful variables
    classLabels = unique(trainY);
    nClasses = length(classLabels);
    params.classLabels = classLabels;
    params.nClasses = nClasses;
    params.nTotal = length(trainY);
    
    
    %% fit LDA model
    for ii=1:nClasses
        params.nClass(ii) = sum(trainY==classLabels(ii));
        params.classPi(ii) = params.nClass(ii) ./ params.nTotal;
        
        curr_class = trainX(trainY==classLabels(ii),:);
        params.Mu(ii,:) = mean(curr_class,1);
        params.Sigma(ii,:,:) = cov(curr_class,1);
    end
    
    params.muTotal = sum(params.Mu,1) ./ params.nClasses;
    params.Sw = squeeze(sum(params.Sigma,1)) ./ params.nClasses;
    params.Sb = zeros(size(params.Sw));
    for ii=1:nClasses
        diffMu = (params.Mu(ii,:) - params.muTotal);
        params.Sb = params.Sb + diffMu'*diffMu;
    end
    params.Sb = params.Sb * params.nTotal;
    
    %% find projection vector and project X 
    invSw_Sb = params.Sw \ params.Sb;
    [u,~,~] = svd(invSw_Sb,'econ');
    
    params.projMat = u;
    params.projVec = u(:,1);
    
    params.projData = trainX*params.projVec;
    
end

