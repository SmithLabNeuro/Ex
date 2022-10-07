function [PD,bl,modVal,modDepth,rSq] = fitCosineTuning(tuningCurves,angleVals)
%FITCOSINETUNING Summary of this function goes here
%   Inputs:
%       - tuningCurves: in Hz (nNeurons x nAngles)
%       - angleVals: angles in degrees (1 x nAngles)
%
%   Outputs:
%       - PD: preferred direction, in degrees (nNeurons x 1)
%       - bl: baseline/average (nNeurons x 1)
%       - modVal: modulation/amplitude of cosine (nNeurons x 1)
%       - modDepth: modulation relative to baseline (nNeurons x 1)
%       - rSq: coefficient of determination of the fit (nNeurons x 1)
%
% Akash Umakantha (aumakant@andrew.cmu.edu)
%

    rad_angleVals = deg2rad(angleVals);
    
    % fit cosine tuning curves
    A = [ones(length(rad_angleVals),1), cos(rad_angleVals), sin(rad_angleVals)];
    betas = A\(tuningCurves');
    
    % get parameters of interest
    PD = atan2(betas(3,:),betas(2,:))';
    bl = betas(1,:)';
    modVal = betas(2,:)'./cos(PD);
    modDepth = modVal./bl;
    
    % assess r-squared
    ss_tot = sum(bsxfun(@minus,tuningCurves,mean(tuningCurves,2)).^2,2);
    tuning_pred = (A*betas)';
    ss_res = sum((tuningCurves-tuning_pred).^2,2);
    rSq = 1 - ss_res./ss_tot;
    
end

