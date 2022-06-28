function [M0, M1, M2, A, Q, C, R, K] = trainKalmanDecoder(binnedSpikesCurrStep, bciValidTrials, joystickKinCurrStep, Qvalue, estParams, latDim, rotMat, separateLatentDimensions, latentKToUse)
  
    allBinnedCountsCurrTime = cat(1, binnedSpikesCurrStep{bciValidTrials})';
    allJoystickKinCurrTime = rotMat' * cat(1, joystickKinCurrStep{bciValidTrials})';
    Tn = cellfun(@(x) size(x, 1), binnedSpikesCurrStep(bciValidTrials));
    Tall = sum(Tn);
    TallLessOne = sum(Tn-1);

    nanTimes = any(isnan(allBinnedCountsCurrTime), 1)...
        | any(isnan(allJoystickKinCurrTime), 1);

    allBinnedCountsCurrTime(:, nanTimes) = [];
    allJoystickKinCurrTime(:, nanTimes) = [];
    Tall = Tall - sum(nanTimes);
    TallLessOne = TallLessOne - sum(nanTimes);

    allLatentProj = latDim * (allBinnedCountsCurrTime - estParams.d);
    meanLatProj = mean(allLatentProj, 2);
    allLatentProj = allLatentProj - meanLatProj;
 
    % compute Kalman parameters for the state model
    A = eye(size(allJoystickKinCurrTime, 1)); %allJoystickKinCurrTime * allJoystickKinPrevTime' / (allJoystickKinPrevTime * allJoystickKinPrevTime');
    Q = Qvalue *  eye(size(allJoystickKinCurrTime, 1)); %[Qvalue 0; 0 Qvalue]; %1/TallLessOne*(allJoystickKinCurrTime - A*allJoystickKinPrevTime) * (allJoystickKinCurrTime - A*allJoystickKinPrevTime)';

    % compute Kalman parameters for the observation model
    C = allLatentProj * allJoystickKinCurrTime' / (allJoystickKinCurrTime * allJoystickKinCurrTime');
    C2 = allLatentProj * allJoystickKinCurrTime(1:2, :)' / (allJoystickKinCurrTime(1:2, :) * allJoystickKinCurrTime(1:2, :)');
%     C3 = allLatentProj * allJoystickKinCurrTime(3, :)' / (allJoystickKinCurrTime(3, :) * allJoystickKinCurrTime(3, :)');
    if separateLatentDimensions
        C = diag(diag(C));
    end
    R = 1/Tall * (allLatentProj - C*allJoystickKinCurrTime) * (allLatentProj - C*allJoystickKinCurrTime)';
    R2 = 1/Tall * (allLatentProj - C2*allJoystickKinCurrTime(1:2, :)) * (allLatentProj - C2*allJoystickKinCurrTime(1:2, :))';
%     R3 = 1/Tall * (allLatentProj - C3*allJoystickKinCurrTime(3, :)) * (allLatentProj - C3*allJoystickKinCurrTime(3, :))';
    if separateLatentDimensions
        R = diag(diag(R));
    end
    % let's try and converge a Kalman gain?
    sigCurrGivenPrev = cov(allJoystickKinCurrTime');
    muCurrGivenPrev = nanmean(allJoystickKinCurrTime, 2);


    for t = 1:100
        Kcurr = sigCurrGivenPrev*C' / (C*sigCurrGivenPrev*C' + R);
        sigCurrGivenCurr = sigCurrGivenPrev - Kcurr*C*sigCurrGivenPrev;
        sigCurrGivenPrev = A*sigCurrGivenCurr*A' + Q;
        Kall{t} = Kcurr;
    end
    
%     sigCurrGivenPrev = cov(allJoystickKinCurrTime(3, :)');
%     muCurrGivenPrev = nanmean(allJoystickKinCurrTime(3, :), 2);
%     
%     
%     for t = 1:100
%         Kcurr = sigCurrGivenPrev*C3' / (C3*sigCurrGivenPrev*C3' + R);
%         sigCurrGivenCurr = sigCurrGivenPrev - Kcurr*C3*sigCurrGivenPrev;
%         sigCurrGivenPrev = A(end, end)*sigCurrGivenCurr*A(end, end)' + Q(end, end);
%         Kall3{t} = Kcurr;
%     end
%     
%     K3 = Kall3{end};
    
    K = Kall{end};
    if separateLatentDimensions && latentKToUse>0
        diagK = diag(K);
        diagK(:) = diagK(latentKToUse);
        K = diag(diagK);
    end

    M1 = A - K*C*A;
    M2 = K * latDim;
    % baseline takes care of accounting for mean offsets in training
    stateModelOffset = zeros(size(A, 1), 1); % for model x_t = Ax_{t-1} + stateModelOffset
    M0 = (A - K*C) * stateModelOffset - M2*estParams.d - K*meanLatProj;
    % M0 = -M1 * joystickKinMean - M2*estParams.d - K*meanLatProj;


end