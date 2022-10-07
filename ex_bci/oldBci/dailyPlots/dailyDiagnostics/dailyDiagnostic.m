function [  ] = dailyDiagnostic( filepath, filename, savePath, savePlots, viewDim )
%DIAGNOSTIC_TRAINING Summary of this function goes here
%   Detailed explanation goes here

    if nargin<4
        savePlots = false;
    end
    if nargin<5
        viewDim = 2;
    end

    % some needed directories
    addpath('bayesClassifier/');
    addpath('../../bciHMM/');
    clc; close all;
    
    % some constants
    targOffCode = 100;
    fixOffCode = 3;
    binSize = 50;
    
    %% decide whether file needs to be loaded or not
    if ~isstruct(filepath)
        load(sprintf('%s/%s',filepath,filename));
    else
        bcistruct = filepath;
    end
    
    %% make directory for saving plots
    tmp = strsplit(filename,'_');
    subjectDate = tmp{1};
    if exist(savePath,'dir')
        warning([savePath ' folder already exists']);
    else
        mkdir(savePath)
    end
    
    %% load parameters and data
    corrTrialCode = 150;
    
    % sorted calibration data
    calibDat = bcistruct.calibrationdata;
    calibDat = calibDat([calibDat.result]==corrTrialCode);
    [~,~,~,calibDat] = prepCalibCounts(calibDat,targOffCode,fixOffCode,0,0.05);
    
    % unsorted calibration data
    unsortedCalibDat = calibDat;
    for ii=1:length(unsortedCalibDat)
        unsortedCalibDat(ii).spikeinfo(:,2) = 0;
    end
    [~,~,~,unsortedCalibDat] = prepCalibCounts(unsortedCalibDat,targOffCode,fixOffCode,0,0.05);
    
    % bci data
    bciDat = bcistruct.bcidata;
    bciDat = getOnlineCounts(bciDat);
    [~,~,~,bciTrialIdx,corrIdx,missIdx] = sepBCIdat(bciDat);
    keepIdx = corrIdx | missIdx;
    bciDat(~keepIdx)=[]; corrIdx(~keepIdx)=[]; missIdx(~keepIdx)=[]; bciTrialIdx(~keepIdx)=[];
    
    % model parameters
    modelparams = bcistruct.modelparamsstruct.modelparams;
    goodNeurons = modelparams.goodneuron;
    %goodNeurons = true(192,1);
    nNeurons = sum(goodNeurons);
    clear('bcistruct');
    
    % some hyperparameters for this day
    bciSessionThresh = bciDat(1).params.block.mindist;
    bciSessionBins = bciDat(1).params.block.updatesOnTarget;
    nDelayBins = modelparams.nDelayBins;
    delayEndIdx = min([calibDat.nBins]);
    nT = delayEndIdx-nDelayBins;
    
    %% get count and target angle data
    % calibration trials
    nCalibTrials = length(calibDat);
    [unsortedCalibCounts,calibCounts] = deal(nan(nNeurons,nT,nCalibTrials));
    for ii = 1:nCalibTrials
        calibCounts(:,:,ii) = calibDat(ii).counts(goodNeurons,(nDelayBins+1):delayEndIdx);
        unsortedCalibCounts(:,:,ii) = unsortedCalibDat(ii).counts(goodNeurons,(nDelayBins+1):delayEndIdx);
    end
    calibAngle = [calibDat.angle];
    uniqueAngles = unique(calibAngle);
    
    % bci trials
    nBciTrials = length(bciDat);
    [unsortedBciCounts,bciCounts] = deal(nan(nNeurons,nT,nBciTrials));
    for ii = 1:nBciTrials
        currT = bciDat(ii).T;
        endIdx = min([currT,nT]);
        unsortedBciCounts(:,1:endIdx,ii) = bciDat(ii).onlineCounts(goodNeurons,1:endIdx);
        bciCounts(:,1:endIdx,ii) = bciDat(ii).onlineSortedCounts(goodNeurons,1:endIdx);
    end
    bciAngle = driftchoiceextractparam(bciDat,'angle=',2);
    uniqueBciAngles = unique(bciAngle);
    
    %% low-d visualization (PCA) of different trial types (different target angles as different colors within each plot)
    % get mean 50 ms spike counts on each trial
    calibTrialMeans = squeeze(nanmean(calibCounts,2))';
    bciTrialMeans = squeeze(nanmean(bciCounts,2))';
    
    % perform pca
    allMeans = [calibTrialMeans; bciTrialMeans];
    [W,pc_score,pc_eigs] = pca(allMeans);
    W = W(:,1:3);
    pc_score = pc_score(:,1:3);
    
    % plot results
    figure; pos=get(gcf,'Position'); set(gcf,'Position',pos.*[1 1 2 2]);
    
    % loadings
    subplot(2,2,1); hold on;
    stem(W(:,1)); stem(W(:,2)); stem(W(:,3));
    xlabel('neuron idx'); ylabel('weight');
    title('PC loadings');
    ylim([-1 1]);
    legend('PC_1','PC_2','PC_3','Location','Best');
    
    % calib
    calib_pcs = pc_score(1:nCalibTrials,:);
    subplot(2,2,2); hold on;
    for i_ang = 1:length(uniqueBciAngles)
        currIdx = calibAngle==uniqueBciAngles(i_ang);
        plot3(calib_pcs(currIdx,1),calib_pcs(currIdx,2),calib_pcs(currIdx,3),'.','MarkerSize',10);
    end
    xlabel('PC_1'); ylabel('PC_2'); zlabel('PC_3');
    legend(num2str(uniqueBciAngles),'Location','Best');
    title('Calibration');
    view(viewDim); xlim([-5 5]); ylim([-5 5]); zlim([-5 5]);
    
    % bci correct
    bci_pcs = pc_score((nCalibTrials+1):end,:);
    subplot(2,2,3); hold on;
    for i_ang = 1:length(uniqueBciAngles)
        currIdx = corrIdx' & bciAngle==uniqueBciAngles(i_ang);
        plot3(bci_pcs(currIdx,1),bci_pcs(currIdx,2),bci_pcs(currIdx,3),'.','MarkerSize',10);
    end
    xlabel('PC_1'); ylabel('PC_2'); zlabel('PC_3');
    legend(num2str(uniqueBciAngles),'Location','Best');
    title('BCI correct');
    view(viewDim); xlim([-5 5]); ylim([-5 5]); zlim([-5 5]);
    
    % bci missed
    subplot(2,2,4); hold on;
    for i_ang = 1:length(uniqueBciAngles)
        currIdx = missIdx' & bciAngle==uniqueBciAngles(i_ang);
        plot3(bci_pcs(currIdx,1),bci_pcs(currIdx,2),bci_pcs(currIdx,3),'.','MarkerSize',10);
    end
    xlabel('PC_1'); ylabel('PC_2'); zlabel('PC_3');
    legend(num2str(uniqueBciAngles),'Location','Best');
    title('BCI missed');
    view(viewDim); xlim([-5 5]); ylim([-5 5]); zlim([-5 5]);
    
    if savePlots
        print([savePath '/pca_trialType.png'],'-dpng');
    end
    
    %% low-d visualization (PCA) of 1st and 2nd trials after a switch vs last trial before a swtich
    % find switches, and other trials of interest around that trial
    firstSwitchIdx = find([0; bciAngle(1:end-1)-bciAngle(2:end)]);
    secondSwitchIdx = min(firstSwitchIdx+1,nBciTrials);
    thirdSwitchIdx = min(firstSwitchIdx+2,nBciTrials);
    
    % plot results
    figure; pos=get(gcf,'Position'); set(gcf,'Position',pos.*[1 1 2 2]);
    
    % first switch trials
    subplot(2,2,1); hold on;
    firstSwitch = false(nBciTrials,1); firstSwitch(firstSwitchIdx) = true;
    firstSwitchAcc = sum(firstSwitch & corrIdx')./sum(firstSwitch);
    for i_ang = 1:length(uniqueBciAngles)
        currIdx = firstSwitch & bciAngle==uniqueBciAngles(i_ang);
        plot3(bci_pcs(currIdx,1),bci_pcs(currIdx,2),bci_pcs(currIdx,3),'.','MarkerSize',10);
    end
    xlabel('PC_1'); ylabel('PC_2'); zlabel('PC_3');
    legend(num2str(uniqueBciAngles),'Location','Best');
    title(sprintf('First switch trial (%.1f %% corr)',firstSwitchAcc*100));
    view(viewDim); xlim([-5 5]); ylim([-5 5]); zlim([-5 5]);
    
    subplot(2,2,2); hold on;
    secondSwitch = false(nBciTrials,1); secondSwitch(secondSwitchIdx) = true;
    secondSwitchAcc = sum(secondSwitch & corrIdx')./sum(secondSwitch);
    for i_ang = 1:length(uniqueBciAngles)
        currIdx = secondSwitch & bciAngle==uniqueBciAngles(i_ang);
        plot3(bci_pcs(currIdx,1),bci_pcs(currIdx,2),bci_pcs(currIdx,3),'.','MarkerSize',10);
    end
    xlabel('PC_1'); ylabel('PC_2'); zlabel('PC_3');
    legend(num2str(uniqueBciAngles),'Location','Best');
    title(sprintf('Second switch trial (%.1f %% corr)',secondSwitchAcc*100));
    view(viewDim); xlim([-5 5]); ylim([-5 5]); zlim([-5 5]);
    
    subplot(2,2,3); hold on;
    thirdSwitch = false(nBciTrials,1); thirdSwitch(thirdSwitchIdx) = true;
    thirdSwitchAcc = sum(thirdSwitch & corrIdx')./sum(thirdSwitch);
    for i_ang = 1:length(uniqueBciAngles)
        currIdx = thirdSwitch & bciAngle==uniqueBciAngles(i_ang);
        plot3(bci_pcs(currIdx,1),bci_pcs(currIdx,2),bci_pcs(currIdx,3),'.','MarkerSize',10);
    end
    xlabel('PC_1'); ylabel('PC_2'); zlabel('PC_3');
    legend(num2str(uniqueBciAngles),'Location','Best');
    title(sprintf('Third switch trial (%.1f %% corr)',thirdSwitchAcc*100));
    view(viewDim); xlim([-5 5]); ylim([-5 5]); zlim([-5 5]);
    
    subplot(2,2,4); hold on;
    nonSwitch = ~(firstSwitch | secondSwitch | thirdSwitch);
    nonSwitchAcc = sum(nonSwitch & corrIdx')./sum(nonSwitch);
    for i_ang = 1:length(uniqueBciAngles)
        currIdx = nonSwitch & bciAngle==uniqueBciAngles(i_ang);
        plot3(bci_pcs(currIdx,1),bci_pcs(currIdx,2),bci_pcs(currIdx,3),'.','MarkerSize',10);
    end
    xlabel('PC_1'); ylabel('PC_2'); zlabel('PC_3');
    legend(num2str(uniqueBciAngles),'Location','Best');
    title(sprintf('Non switch trials (%.1f %% corr)',nonSwitchAcc*100));
    view(viewDim); xlim([-5 5]); ylim([-5 5]); zlim([-5 5]);
    
    if savePlots
        print([savePath '/pca_switchTrials.png'],'-dpng');
    end

%     %% plot tuning curves for each neuron showing each trial and also mean
%     avgCalib = squeeze(mean(calibCounts,2)).*1000/binSize;
%     avgCalib_un = squeeze(mean(unsortedCalibCounts,2)).*1000/binSize;
%     for i_neuron = 1:nNeurons
%         figure; pos=get(gcf,'Position'); set(gcf,'Position',pos.*[0 1 0.5 1]);
%         subplot(2,1,1); hold on; 
%         for i_ang = 1:length(uniqueAngles)
%             currDat = calibTrialMeans(calibAngle==uniqueAngles(i_ang),i_neuron).*1000./binSize;
%             plot(ones(size(currDat)).*uniqueAngles(i_ang),currDat,'r.','MarkerSize',5);
%             plot(uniqueAngles(i_ang),mean(currDat),'r^','MarkerSize',10);
%             xlim([-22.5 360]);
%             xlabel('angle'); ylabel('sp/s');
%             title(sprintf('channel %d calib tuning',i_neuron));
%         end
%         subplot(2,1,2); hold on;
%         plot(1:nCalibTrials,myMovMean(avgCalib(i_neuron,:),[3 3]),'b-');
%         plot(1:nCalibTrials,myMovMean(avgCalib_un(i_neuron,:),[3 3]),'b--');
%         xlabel('calib trial ID'); ylabel('avg trial FR (sp/s)');
%         title(sprintf('channel %d over time',i_neuron));
%         export_fig('~/Desktop/currDat/calibSingleNeuron.pdf','-pdf','-append');
%         close(gcf);
%     end

    %% plot population firing rate (during delay) by trial ID
    % sorted and unsorted calibration means
    popMean_calib = mean(squeeze(mean(calibCounts,2)),1).*1000/binSize;
    popMean_unsortCalib = mean(squeeze(mean(unsortedCalibCounts,2)),1).*1000/binSize;
    
    % sorted and unsorted bci means
    popMean_bci = mean(squeeze(nanmean(bciCounts,2)),1).*1000/binSize;
    popMean_unsortBci = mean(squeeze(nanmean(unsortedBciCounts,2)),1).*1000/binSize;
    
    % smooth data
    nSmoothTrials = 3;
    popMean_calib = myMovMean(popMean_calib,[nSmoothTrials nSmoothTrials]);
    popMean_unsortCalib = myMovMean(popMean_unsortCalib,[nSmoothTrials nSmoothTrials]);
    popMean_bci = myMovMean(popMean_bci,[nSmoothTrials nSmoothTrials]);
    popMean_unsortBci = myMovMean(popMean_unsortBci,[nSmoothTrials nSmoothTrials]);
    
    % plot results
    figure; hold on;
    plot(1:nCalibTrials,popMean_calib,'b-');
    plot(1:nCalibTrials,popMean_unsortCalib,'b--');
    plot((nCalibTrials+1):(nCalibTrials+length(bciDat)),popMean_bci,'r-');
    plot((nCalibTrials+1):(nCalibTrials+length(bciDat)),popMean_unsortBci,'r--');
    xlabel('Trial ID'); ylabel('Pop Avg FR (Hz)'); title('population FR')
    legend('sorted calib','unsorted calib','sorted bci','unsorted bci','Location','Best');
    legend boxoff;
    
    if savePlots
        print([savePath '/popFR.png'],'-dpng');
    end
    
    %% plot distribution of PDs for neurons used on this day
    % compute trial spike counts
    calibTrialCount = squeeze(sum(calibCounts,2));
    
    % compute tuning curves and preferred direction
    [tuneCurves,~,angOrder] = getTuningCurves(calibTrialCount,calibAngle');
    [PDs,bl,modVal,modDepth,rSq] = fitCosineTuning(tuneCurves,angOrder);
    PDs = rad2deg(PDs);
    PDs(PDs<0) = PDs(PDs<0)+360;
    
    % plot results
    figure;
    histogram(PDs,'binWidth',20);
    xlabel('angle'); ylabel('count');
    xlim([0 360]);
    title(sprintf('distribution of PD (%d total neurons)',nNeurons));
    
    if savePlots
        print([savePath '/prefDirDist.png'],'-dpng');
    end
    
    %% confusion matrices during calibration and BCI (bayes decoder)
    % crossvalidate on calibration
    [calib_nb_params,~,calib_pred] = ...
        crossvalidate_nb(calibTrialCount',calibAngle','distType','shareddiag','priorType','equal');
    [calib_confusion,angleOrder] = confusionmat(calib_pred.y,calib_pred.yhat);
    for i_ang = 1:length(angleOrder)
        calib_confusion(i_ang,:) = calib_confusion(i_ang,:) ./ sum(calibAngle==angleOrder(i_ang));
    end
    
    % make plot
    figure;
    imagesc(calib_confusion,[0 .5]); colorbar;
    xlabel('Predicted angle'); ylabel('True angle');
    set(gca,'XTick',find(ismember(angleOrder,uniqueAngles(1:4:end))),'YTick',find(ismember(angleOrder,uniqueAngles(1:4:end))));
    set(gca,'XTickLabel',uniqueAngles(1:4:end),'YTickLabel',uniqueAngles(1:4:end));
    title('Predictions on calib, trained on calib');

    if savePlots
        print([savePath '/calib_confusion.png'],'-dpng');
    end
    
    %% analyze % bci correct over time and as function of target location
    figure; pos=get(gcf,'Position'); set(gcf,'Position',pos.*[1 1 1 2]);

    % percent correct over time
    subplot(2,1,1); hold on;
    bciOutcomes = nan(1,length(bciTrialIdx));
    bciOutcomes(corrIdx==1) = 1;
    bciOutcomes(missIdx==1) = 0;
    trialID = 1:length(bciOutcomes);
    
    for i_ang = 1:length(uniqueBciAngles)
        currIdx = bciAngle==uniqueBciAngles(i_ang);
        tmp = nan(1,length(bciTrialIdx));
        tmp(currIdx) = bciOutcomes(currIdx);
        tmp = myMovMean(tmp,[25 25]);
        plot(trialID,tmp.*100);
    end
    tmp = myMovMean(bciOutcomes,[15 15]);
    plot(trialID,tmp.*100,'k-','LineWidth',3);
    legend(num2str(uniqueBciAngles),'Location','Best'); legend boxoff;
    xlim([1 length(trialID)]); ylim([0 100]);
    xlabel('BCI trial ID'); ylabel('BCI % correct');
    title('BCI performance over time (smoothed over 30 trials)');
    
    % percent correct as a function of target direction
    subplot(2,1,2); hold on;
    angCorr = nan(size(uniqueBciAngles));
    for i_ang = 1:length(uniqueBciAngles)
        angCorrIdx = bciOutcomes(bciAngle==uniqueBciAngles(i_ang));
        angCorr(i_ang) = sum(angCorrIdx)./length(angCorrIdx);
    end
    plot(uniqueBciAngles,angCorr.*100,'o-');
    xlabel('Target Location'); ylabel('BCI % correct');
    xlim([min(uniqueBciAngles)-10 max(uniqueBciAngles)+10]); ylim([-1 101]);
    set(gca,'XTick',uniqueBciAngles);
    
    if savePlots
        print([savePath '/bciPerformance.png'],'-dpng');
    end
    
    
    %% analysis of different reward criterion on bci data
    
%     if savePlots
%         print([mySavePath '/rewardCriterion.png'],'-dpng');
%     end
    
end

