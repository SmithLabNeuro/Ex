function [ ] = stabilityDiagnostic(filepath, filename, savePath, savePlots)
%STABILITYDIAGNOSTIC Summary of this function goes here
%   Detailed explanation goes here

    close all;

    if nargin<4
        savePlots = false;
    end
    
    % some constants
    fixOnCode = 140;
    fixOffCode = 3;
    binSize = 50;
    bciCorrCode = 161;
    bciMissCode = 162;
    corrTrialCode = 150;
    
    % some parameters
    nSmooth = 25;
    
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
    saveName = sprintf('%s/stabilityDiagnostic_%s.pdf',savePath,subjectDate);
    
    %% load parameters and data
    % get model params
    modelparams = bcistruct.modelparamsstruct(end).modelparams;
    faParams = modelparams.estParams;
    zDim = size(faParams.L,2);
    channelIdx = 1:192;
    channelIdx = channelIdx(modelparams.goodneuron);
    nNeurons = length(channelIdx);
    
    % sorted cross-talk data
    crosstalkDat = bcistruct.calibrationdata;
    crosstalkDat = crosstalkDat([crosstalkDat.result]==corrTrialCode);
    [~,~,~,crosstalkDat] = prepCalibCounts(crosstalkDat,fixOnCode,fixOffCode,0,binSize/1000);
    
    % get fa calibration data, and bci data
    bciDat = bcistruct.bcidata;
    bciDat = getOnlineCounts(bciDat);
    bciDat = bciDat(ismember([bciDat.result],[bciCorrCode bciMissCode]));
    calibDat = bciDat(1:60);
    bciDat = bciDat(61:end);
    nT = min([calibDat.T]);
    t = (1:nT).*50 - 50;
    
    %% get counts
    % cross-talk trials
    nXTalkTrials = length(crosstalkDat);
    xTalkCounts = deal(nan(nNeurons,nT,nXTalkTrials));
    xTalkIdx = 1:nXTalkTrials;
    for ii=1:nXTalkTrials
        xTalkCounts(:,:,ii) = crosstalkDat(ii).counts(modelparams.goodneuron,1:nT);
    end
    
    % calibration trials
    nCalibTrials = length(calibDat);
    calibCounts = deal(nan(nNeurons,nT,nCalibTrials));
    calibIdx = (nXTalkTrials+1):(nXTalkTrials+nCalibTrials);
    for ii = 1:nCalibTrials
         tmp = calibDat(ii).onlineSortedCounts;
         currBins = min(size(tmp,2),nT);
         calibCounts(:,1:currBins,ii) = tmp(modelparams.goodneuron(1:size(tmp,1)),1:currBins);
    end
    
    % bci trials
    nBciTrials = length(bciDat);
    bciCounts = deal(nan(nNeurons,nT,nBciTrials));
    bciIdx = (nXTalkTrials+nCalibTrials+1):(nXTalkTrials+nCalibTrials+nBciTrials);
    for ii = 1:nBciTrials
         tmp = bciDat(ii).onlineSortedCounts;
         currBins = min(size(tmp,2),nT);
         bciCounts(:,1:currBins,ii) = tmp(modelparams.goodneuron(1:size(tmp,1)),1:currBins);
    end
    
    %% percent shared of each neuron
    % actual
    L = modelparams.estParams.L; Ph = modelparams.estParams.Ph;
    sharedVar = diag(L*L');
    percShared = sharedVar./(sharedVar+Ph) .* 100;
    
    figure;
    stem(channelIdx,percShared,'b');
    xlabel('channel #'); ylabel('% shared var');
    
    if savePlots
        export_fig(saveName,'-pdf','-painters','-append','-nocrop');
    end
    
    %% percent correct and acquisition time
    calibCorr = [calibDat.result]==bciCorrCode;
    calibAcqTime = [calibDat.T];
    calibAcqTime(~calibCorr) = nan;
    
    bciCorr = [bciDat.result]==bciCorrCode;
    bciAcqTime = [bciDat.T];
    bciAcqTime(~bciCorr) = nan;
    
    figure;
    
    yyaxis left; hold on;
    plot(calibIdx,myMovMean(calibAcqTime,[nSmooth nSmooth]).*binSize,'b-');
    plot(bciIdx,myMovMean(bciAcqTime,[nSmooth nSmooth]).*binSize,'b-');
    xlabel('trial #'); ylabel('acq time (ms)');
    yyaxis right; hold on;
    plot(calibIdx,myMovMean(calibCorr,[nSmooth nSmooth]).*100,'r-');
    plot(bciIdx,myMovMean(bciCorr,[nSmooth nSmooth]).*100,'r-');
    ylim([0 100]);
    ylabel('% correct');
    plot([1 1].*nXTalkTrials,ylim,'k--','LineWidth',2);
    plot([1 1].*(nXTalkTrials+nCalibTrials),ylim,'k--','LineWidth',2);
    xlim([1 bciIdx(end)]);
    title(sprintf('Behavior (smoothed over %d trials)',2*nSmooth));
    
    if savePlots
        export_fig(saveName,'-pdf','-painters','-append','-nocrop');
    end
    
    %% population firing rate over course of the session (i.e. neural drifts)
    popMean_xTalk = nanmean(squeeze(nanmean(xTalkCounts,1)),1).*1000/binSize;
    popMean_calib = nanmean(squeeze(nanmean(calibCounts,1)),1).*1000/binSize;
    popMean_bci = nanmean(squeeze(nanmean(bciCounts,1)),1).*1000/binSize;
    
    figure; hold on;
    plot(xTalkIdx,myMovMean(popMean_xTalk,[nSmooth nSmooth]),'g-');
    plot(calibIdx,myMovMean(popMean_calib,[nSmooth nSmooth]),'b-');
    plot(bciIdx,myMovMean(popMean_bci,[nSmooth nSmooth]),'r-');
    plot([1 1].*nXTalkTrials,ylim,'k--','LineWidth',2);
    plot([1 1].*(nXTalkTrials+nCalibTrials),ylim,'k--','LineWidth',2);
    xlabel('trial #'); ylabel('avg pop firing rate (Hz)');
    title(sprintf('Population rate (smoothed over %d trials)',2*nSmooth));
    
    if savePlots
        export_fig(saveName,'-pdf','-painters','-append','-nocrop');
    end
    
    %% plot distribution of latent projects before and after smoothing
    figure; pos=get(gcf,'Position'); set(gcf,'Position',pos.*[1 1 3 2]);
    
    xTalkFAproj = nan(zDim,size(xTalkCounts,2),nXTalkTrials);
    for i_trial = 1:nXTalkTrials
        tmp = fastfa_estep(xTalkCounts(:,:,i_trial),faParams);
        xTalkFAproj(:,:,i_trial) = expSmoothing(tmp.mean,modelparams.alpha);
    end
    
    calibFAproj = nan(zDim,size(calibCounts,2),nCalibTrials);
    for i_trial = 1:nCalibTrials
        tmp = fastfa_estep(calibCounts(:,:,i_trial),faParams);
        calibFAproj(:,:,i_trial) = expSmoothing(tmp.mean,modelparams.alpha);
    end
    
    for i_dim = 1:zDim
        subplot(2,zDim,i_dim); hold on;
        tmp = squeeze(calibFAproj(i_dim,:,:));
        histogram(tmp(:),'BinWidth',0.1,'Normalization','probability'); 
        plot([1 1].*modelparams.meanlatent(i_dim),ylim,'k--','LineWidth',2);
        xlabel('proj value'); xlim([-2.5 2.5])
        title(sprintf('smoothed FA_%d, calib',i_dim));
        legend('','\mu_{calib}','Location','Best');
    end
    
    bciFAproj = nan(zDim,size(bciCounts,2),nBciTrials);
    bciFAproj_unsmooth = nan(zDim,size(bciCounts,2),nBciTrials);
    for i_trial = 1:nBciTrials
        tmp = fastfa_estep(bciCounts(:,:,i_trial),faParams);
        bciFAproj_unsmooth(:,:,i_trial) = tmp.mean;
        bciFAproj(:,:,i_trial) = expSmoothing(tmp.mean,modelparams.alpha);
    end
    
    for i_dim = 1:zDim
        subplot(2,zDim,i_dim+zDim); hold on;
        tmp = squeeze(bciFAproj(i_dim,:,:));
        histogram(tmp(:),'BinWidth',0.1,'Normalization','probability'); 
        plot([1 1].*modelparams.meanlatent(i_dim),ylim,'k--','LineWidth',2);
        plot([1 1].*nanmean(tmp(:)),ylim,'r--','LineWidth',2);
        xlabel('proj value'); xlim([-2.5 2.5])
        title(sprintf('smoothed FA_%d, bci',i_dim));
        legend('','\mu_{calib}','\mu_{bci}','Location','Best');
    end
    
    if savePlots
        export_fig(saveName,'-pdf','-painters','-append','-nocrop');
    end
    
    %% latents over course of the session (i.e. latent drifts)
    figure; pos=get(gcf,'Position'); set(gcf,'Position',pos.*[1 1 3 1]);
    
    xTalkFAavg = squeeze(nanmean(xTalkFAproj,2));
    calibFAavg = squeeze(nanmean(calibFAproj,2));
    bciFAavg = squeeze(nanmean(bciFAproj,2));
    
    for i_dim = 1:zDim
        subplot(1,zDim,i_dim); hold on;
        plot(xTalkIdx,myMovMean(xTalkFAavg(i_dim,:),[nSmooth nSmooth]),'g-');
        plot(calibIdx,myMovMean(calibFAavg(i_dim,:),[nSmooth nSmooth]),'b-');
        plot(bciIdx,myMovMean(bciFAavg(i_dim,:),[nSmooth nSmooth]),'r-');
        ylim([-0.4 0.4]);
        plot([1 1].*nXTalkTrials,ylim,'k--','LineWidth',2);
        plot([1 1].*(nXTalkTrials+nCalibTrials),ylim,'k--','LineWidth',2);
        plot(xlim,[1 1].*modelparams.meanlatent(i_dim),'k--','LineWidth',2);
        xlabel('trial #'); ylabel('avg proj');
        title(sprintf('Avg FA_%d (smoothed over %d trials)',i_dim,2*nSmooth));
    end
    
    if savePlots
        export_fig(saveName,'-pdf','-painters','-append','-nocrop');
    end
    
    %% plot distribution of distances from mean & chose threshold value (calib and bci)
    figure; hold on;
    
    xTalkDist = squeeze(sum(bsxfun(@minus,xTalkFAproj,modelparams.meanlatent).^2,1));
    calibDist = squeeze(sum(bsxfun(@minus,calibFAproj,modelparams.meanlatent).^2,1));
    bciDist = squeeze(sum(bsxfun(@minus,bciFAproj,modelparams.meanlatent).^2,1));
    
    calibDistHist = calibDist(11:end,:); calibDistHist = calibDistHist(:);
    bciDistHist = bciDist(11:end,:); bciDistHist = bciDistHist(~isnan(bciDistHist));
    histogram(calibDistHist,'BinWidth',0.1,'Normalization','probability','DisplayStyle','stairs');
    histogram(bciDistHist,'BinWidth',0.1,'Normalization','probability','DisplayStyle','stairs');
    plot([1 1].*modelparams.threshold,ylim,'k--');
    xlabel('sq distance'); title('sq distances');
    legend('calib','bci','distance thresh','Location','Best');
    xlim([0 4]);
    
    if savePlots
        export_fig(saveName,'-pdf','-painters','-append','-nocrop');
    end
    
    %% plot distance on each trial over course of the session (distance drifts)
    figure; hold on;
    
    xTalkDistAvg = nanmean(xTalkDist(11:end,:),1);
    calibDistAvg = nanmean(calibDist(11:end,:),1);
    bciDistAvg = nanmean(bciDist(11:end,:),1);
    
    plot(xTalkIdx,myMovMean(xTalkDistAvg,[nSmooth nSmooth]),'g-');
    plot(calibIdx,myMovMean(calibDistAvg,[nSmooth nSmooth]),'b-');
    plot(bciIdx,myMovMean(bciDistAvg,[nSmooth nSmooth]),'r-');
    plot(xlim,[1 1].*modelparams.threshold,'k--','LineWidth',2);
    plot([1 1].*nXTalkTrials,ylim,'k--','LineWidth',2);
    plot([1 1].*(nXTalkTrials+nCalibTrials),ylim,'k--','LineWidth',2);
    xlabel('trial #'); ylabel('avg sq dist');
    title(sprintf('Avg sq dist (smoothed over %d trials)',2*nSmooth));
    
    if savePlots
        export_fig(saveName,'-pdf','-painters','-append','-nocrop');
    end
    
    %% plot smoothed hit rates vs smoothed FA projections
    figure; pos=get(gcf,'Position'); set(gcf,'Position',pos.*[1 1 3 1.5]);
    
    smoothHit = myMovMean(bciCorr,[nSmooth nSmooth]).*100;
    bciFAavg = squeeze(nanmean(bciFAproj(:,21:end,:),2));
    for i_dim = 1:zDim
        tmp = myMovMean(bciFAavg(i_dim,:),[nSmooth nSmooth]);
        subplot(2,zDim,i_dim); hold on;
        plot(smoothHit,tmp,'-');
        xlabel('hit rate'); ylabel(sprintf('Avg FA_%d',i_dim));
        title(sprintf('smoothed over %d trials, corr=%.2f',2*nSmooth,corr(smoothHit',tmp')));
        xlim([20 100]); ylim([-.3 .3]);
    end
    
    bciDist_dim = bsxfun(@minus,bciFAproj,modelparams.meanlatent).^2;
    bciDist_dim = squeeze(nanmean(bciDist_dim(:,21:end,:),2));
    for i_dim = 1:zDim
        tmp = myMovMean(bciDist_dim(i_dim,:),[nSmooth nSmooth]);
        subplot(2,zDim,i_dim+zDim); hold on;
        plot(smoothHit,tmp,'-');
        xlabel('hit rate'); ylabel(sprintf('Avg Sq Dist FA_%d',i_dim));
        title(sprintf('smoothed over %d trials, corr=%.2f',2*nSmooth,corr(smoothHit',tmp')));
    end
    
    if savePlots
        export_fig(saveName,'-pdf','-painters','-append','-nocrop');
    end
    
    %% plot single trial trajectories on miss trials and correct trials
    missIdx = [bciDat.result]==162;
    bciMissDist = bciDist(:,missIdx);
    bciCorrDist = bciDist(:,~missIdx);
    
    
    figure; hold on;
    plot(t,nanmean(bciMissDist,2),'b-');
    plot(xlim,[1 1].*modelparams.threshold,'k--','LineWidth',2);
    xlabel('time (ms)')
    ylabel('euc dist');
    title('Miss trials');
    xlim([500 3500]);
    
    %% plot miss trial trajectories along each FA dimension
    missTrials = bciFAproj(:,:,missIdx);
    missTrials_unsmooth = bciFAproj_unsmooth(:,:,missIdx);
    
    figure; pos=get(gcf,'Position'); set(gcf,'Position',pos.*[1 1 2 1]);
    for i_dim = 1:zDim
        tmp_unsmooth = squeeze(missTrials_unsmooth(i_dim,:,:));
        tmp = squeeze(missTrials(i_dim,:,:));
        subplot(1,2,1); hold on;
        plot(t,nanmean(tmp,2),'-');
        subplot(1,2,2); hold on;
        plot(t,nanmean(tmp_unsmooth,2),'-');    
    end
    subplot(1,2,1); hold on;
    xlabel('time (ms)'); ylabel('avg proj');
    title('smoothed FA proj (miss trials)'); xlim([500 3500]);
    legend('FA_1','FA_2','FA_3','FA_4','FA_5','Location','Best');
    subplot(1,2,2); hold on;
    xlabel('time (ms)'); ylabel('avg proj');
    title('raw FA proj (miss trials)'); xlim([500 3500]);
    legend('FA_1','FA_2','FA_3','FA_4','FA_5','Location','Best');
    
    figure; pos=get(gcf,'Position'); set(gcf,'Position',pos.*[1 1 3 1]);
    rngIdx = randperm(sum(missIdx));
    rngIdx = rngIdx(1:20);
    for i_dim = 1:zDim
        subplot(1,zDim,i_dim); hold on;
        tmp = squeeze(missTrials(i_dim,:,:));
        plot(t,tmp(:,rngIdx),'b-');
        xlabel('trial #'); ylabel('avg proj');
        title(sprintf('FA_%d',i_dim));
    end
    
end

