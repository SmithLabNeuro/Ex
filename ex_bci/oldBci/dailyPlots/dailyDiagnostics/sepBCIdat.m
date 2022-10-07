function [ xMGS, xBciCorr, xBciMiss, bciTrialIdx, corrIdx, missIdx ] = sepBCIdat( bciDat )
%SEPBCIDAT Summary of this function goes here
%   Detailed explanation goes here

    nTrials = length(bciDat);
    
    
    bciCorrectCode = 161;
    bciMissCode = 162;
    
    for i_trial = 1:nTrials
        bciDat(i_trial).bciTrial = bciDat(i_trial).params.trial.bciTrial;
        bciDat(i_trial).bciCorr = sum(bciDat(i_trial).trialcodes(:,2)==bciCorrectCode)>0;
        bciDat(i_trial).bciMiss = sum(bciDat(i_trial).trialcodes(:,2)==bciMissCode)>0;
    end

    % separate into the proper structures
    bciTrialIdx = [bciDat.bciTrial];
    corrIdx = [bciDat.bciCorr];
    missIdx = [bciDat.bciMiss];
    xMGS = bciDat(bciTrialIdx==0);
    xMGS = xMGS([xMGS.result]==150);
    xBciCorr = bciDat(corrIdx==1);
    xBciMiss = bciDat(missIdx==1);

    fprintf('Percent of BCI correct trials: %.2f\n',...
        sum(corrIdx==1)/sum(corrIdx==1 | missIdx==1)*100);
    
end

