function [ params, trainError ] = update_nb( params, trainX, trainY )
%UPDATE_NB Summary of this function goes here
%   Detailed explanation goes here

    % update sufficient statistics
    for ii = 1:params.nClasses
        curr_class = trainX(trainY==params.classLabels(ii),:);
        params.suffStat(ii).sumX = params.suffStat(ii).sumX + sum(curr_class,1);
        params.suffStat(ii).N = params.suffStat(ii).N + size(curr_class,1);
        params.suffStat(ii).sumX2 = params.suffStat(ii).sumX2 + sum(curr_class.^2,1);
    end
    
    % update class priors
    switch lower(params.priorType)
        case 'setprior'
            for ii = 1:params.nClasses
                params.classPi(ii) = params.suffStat(ii).N./sum([params.suffStat.N]);
            end
        case 'equalprior'
            for ii = 1:params.nClasses
                params.classPi(ii) = 1./params.nClasses;
            end
        otherwise
            error('Please specify a valid prior type')
    end  
    
    
    % update model parameters
    switch lower(params.distType)
        case 'diaggauss' % each class gets its own cov, but it's diagonal
            % fit means
            for ii=1:params.nClasses
                params.Mu(ii,:) = params.suffStat(ii).sumX./params.suffStat(ii).N;
            end
            % fit covariance
            for ii=1:params.nClasses
                E_X2 = params.suffStat(ii).sumX2./(params.suffStat(ii).N);
                E_X_2 = params.Mu(ii,:).^2;
                covX = E_X2 - E_X_2;
                params.Sigma(ii,:,:) = diag(max(covX,params.minVar));
            end
        case 'fullgauss' % each class gets its own cov
            error('Cannot update fullgauss model')
        case 'sharedgauss' % gaussian w/ same diagonal cov for each class
            % fit means
            for ii=1:params.nClasses
                params.Mu(ii,:) = params.suffStat(ii).sumX./params.suffStat(ii).N;
            end
            % fit covariance
            modelCov = zeros(size(trainX,2));
            for ii=1:params.nClasses
                E_X2 = params.suffStat(ii).sumX2./(params.suffStat(ii).N);
                E_X_2 = params.Mu(ii,:).^2;
                covX = E_X2 - E_X_2;
                modelCov = modelCov + params.classPi(ii) .* diag(covX);
            end
            for ii=1:params.nClasses
                % diagonalize the matrix and set a minimum variance
                params.Sigma(ii,:,:) = diag(max(diag(modelCov),params.minVar));
            end
        case 'poisson'
            % fit lambda
            for ii=1:params.nClasses
                params.Lambda(ii,:) = params.suffStat(ii).sumX./params.suffStat(ii).N;
            end
            % set a minimum lambda
            params.Lambda = max(params.Lambda,params.minVar);
        otherwise
            error('Please specify a valid distribution type')
    end

    % calculate training error
    predictions = predict_nb(trainX,params);
    trainError = sum(predictions~=trainY)./length(trainY);
    
end

