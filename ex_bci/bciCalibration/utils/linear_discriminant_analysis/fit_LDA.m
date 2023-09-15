function [params] = fit_LDA(trainX, trainY)
% Fits a linear discriminant model to provided trainX and trainY
% Based off of code written by Akash Umakantha: https://github.com/akash-uma/multiDim_lda/blob/master/train_lda.m

    % useful variables
    classLabels = unique(trainY);
    nClasses = length(classLabels);
    params.classLabels = classLabels;
    params.nClasses = nClasses;
    params.nTotal = length(trainY);
    
    
    %% fit LDA model
    for ii=1:nClasses
        % Keep track of class member totals
        params.nClass(ii) = sum(trainY==classLabels(ii));
        % Keep track of proportion of each class members to entire
        % population
        params.classPi(ii) = params.nClass(ii) ./ params.nTotal;
        
        curr_class_X = trainX(trainY==classLabels(ii),:);
        % Keep track of class means and covariances 
        params.Mu(ii,:) = mean(curr_class_X,1);
        % Normalizes by N instead of N-1 (Bessel's correction not applied
        % here). 
        params.Sigma(ii,:,:) = (curr_class_X-params.Mu(ii,:))'*(curr_class_X-params.Mu(ii,:));
    end
    
    % Take mean of total training set
    params.muTotal = mean(trainX,1);
    % Generate within-scatter matrix (should be the sum of product of
    % covariance matrices and their respective num_class_members)
    params.Sw = squeeze(sum(params.Sigma,1)) ;
    params.Sb = zeros(size(params.Sw));
    % Generate between-class scatter matrix 
    for ii=1:nClasses
        diffMu = (params.Mu(ii,:) - params.muTotal);
        % Sum up the outer product of these differences
        params.Sb = params.Sb + params.nClass(ii)*(diffMu'*diffMu);
    end
    %% find projection vector and project X 
    % identify solution to A*V = B*V*D where V is our LDA axes and D are
    % the eigenvalues.
    [V,D] =  eig(params.Sb, params.Sw);
    params.projMat = V;
    params.projVec = V(:,1);
    params.projData = trainX*params.projVec;
end