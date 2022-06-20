function modelParams = train_lowd_kf(trainDat,modelParams,varargin)

    %% set optional arguments
    nDelayBins = 4;
    smooth_scale = 5e-3;
    prior_strength = 0;
    assignopts(who,varargin);
    
    %% get some important variables
    nNeurons = size(trainDat(1).counts,1);
    nTrials = length(trainDat);
    minBins = min([trainDat.nBins]);
    nTimePoints = minBins - nDelayBins;
    targAngles = [trainDat.angle]';
    modelParams.allangles = unique(targAngles);
    
    %% parse trainDat into a matrix of nNeurons x nTrials x nTimePoints
    allCounts = nan(nNeurons,nTimePoints,nTrials);
    startIdx = nDelayBins+1;
    endIdx = nDelayBins+nTimePoints;
    for i_trial = 1:nTrials
        allCounts(:,:,i_trial) = trainDat(i_trial).counts(:,startIdx:endIdx);
    end
    
    % only keep good neurons
    allCounts = allCounts(modelParams.goodneuron,:,:);
    nNeurons = size(allCounts,1);
    
    %% project into FA space
    lowd_x = nan(modelParams.zDim,nTimePoints,nTrials);
    for i_trial = 1:nTrials
        tmp = fastfa_estep(squeeze(allCounts(:,:,i_trial)),modelParams.fa_params);
arnew.sig = sig_1_0 - kt*C*sig_1_        tmp = orthogonalize(tmp.mean,modelParams.fa_params.L);
        lowd_x(:,:,i_trial) = tmp;
    end
    
    %% generate cursor positions
    targ_comp = exp(deg2rad(modelParams.allangles)*1i);
    targ_xy = [real(targ_comp), imag(targ_comp)]';
    
    train_pos = nan(2,size(lowd_x,2),nTrials);
    
    for ii = 1:nTrials
        curr_xy = targ_xy(:,targAngles(ii)==modelParams.allangles);
        train_pos(1,:,ii) = curr_xy(1);
        train_pos(2,:,ii) = curr_xy(2);
    end
    
    %% train kalman filter
    tmp = reshape(lowd_x,size(lowd_x,1),[]);
    kf_params.d = mean(tmp,2);
    lowd_x = bsxfun(@minus,lowd_x,kf_params.d);
    
    zt = reshape(train_pos,size(train_pos,1),[]);
    xt = reshape(lowd_x,size(lowd_x,1),[]);
    
    kf_params.pival = [0; 0];
    kf_params.V = [1 0; 0 1].*prior_strength;
    
    kf_params.A = [1 0; 0 1];
    kf_params.Q = [1 0; 0 1].*smooth_scale;
    
    C = (xt*zt')*inv(zt*zt');
    resid_vals = xt - C*zt;
    
    kf_params.C = C;
    kf_params.R = (resid_vals*resid_vals')./size(resid_vals,2);
    
    modelParams.kf_params = kf_params;
    
end

