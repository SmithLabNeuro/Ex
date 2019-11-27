function [modelparams,modelparams1,modelparams2] = calibrateDistanceBCIFA(filepath,savepath,varargin)
% [modelparams, outdat] = calibrateDistanceBCIFA('D:/pepelepew/','X:/pepelepew/','gamma',0.2,'targetcorrect',0.5,'zdim',4)
addpath bcihelpers
filename=[];
gamma =0.2; % gamma value for neural network sorting
targetCorrect = 0.5; % the target baseline correct for the bci system
offlinemodeflag=0;
trialnum=0; 
dat=[];
modelparamsold=[];
calibbin=0.05; % size of count bin
numbins=8; % number of bins needed for bci correct
waitbins=8; % number of bins to wait after bci start
zdim=4; % number of fa dimensions
newcalcflag=1; % this is typically set to one during online calibration
dimsbci1 = [1 2 3 4];
dimsbci2 = [];
highscale = 5;
deltabciflag = 0;
assignopts(who,varargin);
display(['Using deltabciflag = ', num2str(deltabciflag)])

%filename = [];
if offlinemodeflag == 0
%% find file
    if trialnum==0
        if isempty(filename)
            filename = findBCIFileName(filepath);
        end

        %% make folder for the day
        tmp = strsplit(filename,'_');
        subjectDate = tmp{1};
        %mySavePath = ['C:/Users/smithlab/dropbox/smithlabdata/bcidailyplots/' subjectDate];
        mySavePath = [savepath subjectDate];
        if exist(mySavePath,'dir')
            warning([mySavePath ' folder already exists'])
        else
            mkdir(mySavePath)
        end
        addpath(genpath('../structBuilders/'));
        dat = nev2sortedStruct([filepath,filename,'.nev'],gamma,0);
    else
        trainingdat = dat;
    end
else

    %% offline loading
    filename = [filepath,filename];%'../Pe180906bcistruct';
    load(filename)
    addpath('C:/Users/rwill/Dropbox/smithlabrig/Ex/ex_bci/calib_online/');
end


if trialnum == 0||newcalcflag ==0
    countbinsize = calibbin;
    modelparams.targetCorrect = targetCorrect;
    modelparams.countbinsize = countbinsize;
    modelparams.numbins = numbins;
    modelparams.waitbins = waitbins;
    waitbins = waitbins - 1; % 2018/11/28, RW do this b/c there is a communication delay between control and bci computer so the bci task only waits for waitbins-1 bins
    modelparams.highscale = highscale;
    results = [dat.result];
    trainingdat = dat(results==150);
    onlythesetrials = 1:length(trainingdat);
    modelparams.onlythesetrials = onlythesetrials;
    trainingdat = trainingdat(onlythesetrials);
    [countmatpre,trainingdat] = prepCalibCounts(trainingdat,140, 3,countbinsize);
    
    %% determine good neurons
    meanRate = mean(countmatpre,2).*(1/countbinsize);
    fanoFact = var(countmatpre,[],2)./mean(countmatpre,2);
    goodneuron = meanRate>1 & (fanoFact<3 & fanoFact>0.5);%meanRate > 2 & fanoFact<3 & fanoFact>0.5;
    goodneuron(97:end) = false;
    channelIdx = 1:192;
    channelIdx = channelIdx(goodneuron);
    nNeurons = sum(goodneuron);

    %% find cross-talk
    % get 1 ms counts
    [~,tmpDat] = prepCalibCounts(trainingdat,140,3,0.001);
    minBins = floor(min([tmpDat.nBins])/1000)*1000;
    spTrains = nan(sum(goodneuron),minBins,length(tmpDat));
    for ii=1:length(tmpDat)
        spTrains(:,:,ii) = tmpDat(ii).counts(goodneuron,1:minBins);
    end
    
    % shuffle control
    spTrainShuffle = nan(size(spTrains)); 
    for ii=1:nNeurons
        rngIdx = randperm(size(spTrains,3));
        spTrainShuffle(ii,:,:) = spTrains(ii,:,rngIdx);
    end
    nTrials = size(spTrains,3);
    
    % compute trial ccgs
    nLag = 150;
    modelparams.nLag = nLag;
    xCorr_trial = nan(nNeurons,nNeurons,nLag*2+1,nTrials);
    xCorr_shuffleTrial = nan(size(xCorr_trial));
    for i_trial = 1:nTrials
        fprintf('      CCGs on trial %d of %d...\n',i_trial,nTrials);
        currDat = spTrains(:,:,i_trial);
        xCorr_trial(:,:,:,i_trial) = myXcorr(currDat,nLag,false);
        shuffleDat = spTrainShuffle(:,:,i_trial);
        xCorr_shuffleTrial(:,:,:,i_trial) = myXcorr(shuffleDat,nLag,false);
    end
    avg_xcorr = squeeze(mean(xCorr_trial,4));
    avg_shuffleXcorr = squeeze(mean(xCorr_shuffleTrial,4));
    diff_xcorr = avg_xcorr - avg_shuffleXcorr;
    
    % find large lag 0 correlations
    statMat = diff_xcorr; 
    statMat(:,:,nLag:(nLag+2)) = [];
    noLagVals = diff_xcorr(:,:,nLag+1);
    meanVals = mean(statMat,3);
    stdVals = std(statMat,1,3);
    threshVals = meanVals+6*stdVals;
    modelparams.threshVals = threshVals;
    % remove channels with crosstalk
    [x,y] = find(noLagVals>threshVals);
    while ~isempty(x)
        rmNeuron = mode([x; y]);
        noLagVals(rmNeuron,:) = []; noLagVals(:,rmNeuron) = [];
        threshVals(rmNeuron,:) = []; threshVals(:,rmNeuron) = [];
        channelIdx(rmNeuron) = [];

        [x,y] = find(noLagVals>threshVals);
    end

%     % get spike train matrix
%     [~,~,~,tmpDat,~,~] = prepCalibCounts(trainingdat,140,3,0,0.0001,[]);
%     minBins = floor(min([tmpDat.nBins])/1000)*1000;
%     spTrains = nan(sum(goodneuron),minBins,length(tmpDat));
%     for ii=1:length(tmpDat)
%         spTrains(:,:,ii) = tmpDat(ii).counts(goodneuron,1:minBins);
%     end
%     fullSpTrain = reshape(spTrains,sum(goodneuron),[]);
%     
%     % compute normalized coincidences
%     c = normCoincidence(fullSpTrain);
%     
%     % remove channels with crosstalk
%     [x,y] = find(c>=0.1);
%     while ~isempty(x)
%         rmNeuron = mode([x; y]);
%         c(rmNeuron,:) = []; c(:,rmNeuron) = [];
%         channelIdx(rmNeuron) = [];
% 
%         [x,y] = find(c>=0.1);
%     end

    goodneuron = false(192,1);
    goodneuron(channelIdx) = true;
    
    fprintf('   %d of %d are good neurons\n',sum(goodneuron),96);
    
    countmat = countmatpre(goodneuron,:);
    modelparams.goodneuron = goodneuron;
    modelparams.gamma = gamma;
    for n = 1:length(trainingdat)
        trainingdat(n).counts = trainingdat(n).counts(modelparams.goodneuron,:);
    end
else
    countbinsize = modelparamsold(1).countbinsize;
    modelparams.targetCorrect = modelparamsold(1).targetCorrect;
    modelparams.countbinsize = modelparamsold(1).countbinsize;
    modelparams.numbins = modelparamsold(1).numbins;
    modelparams.waitbins = modelparamsold(1).waitbins;
    waitbins = waitbins - 1; % 2018/11/28, RW do this b/c there is a communication delay between control and bci computer so the bci task only waits for waitbins-1 bins
    modelparams.highscale = modelparamsold(1).highscale;
    
    
    countmat = [];
    modelparams.goodneuron = modelparamsold(1).goodneuron;
    modelparams.gamma = modelparamsold(1).gamma;
    goodneuron = modelparams.goodneuron;%ones(size(countmat,1),1)==1;
    for n = 1:length(trainingdat)
        countmat = [countmat trainingdat(n).counts];
    end
    
    %the following is needed to prevent neurons dropped during calibration
    %iterations from causing a low-rank covariance matrix, which would
    %crash the fastfa.m function.

%     meanRate = mean(countmat,2).*(1/countbinsize);
%     fanoFact = var(countmat,[],2)./mean(countmat,2);
%     thisinds = goodneuron==1;
%     goodneuron(goodneuron==1) = meanRate > 2 & fanoFact<2 & fanoFact>0.5;
%     countmat = countmat(goodneuron(thisinds),:);
%     for n = 1:length(trainingdat)
%         trainingdat(n).counts = trainingdat(n).counts(goodneuron(thisinds),:);
%     end
%     modelparams.goodneuron = goodneuron;
end

%% define mu %%% MODEL SPECIFIC CODE %%%
modelparams.filtertimeconstanst = 6;
modelparams.alpha = 1-exp(-1/modelparams.filtertimeconstanst);
modelparams.mu = mean(countmat,2);
if deltabciflag == 1
    distancemetric = @(x1,x2,mu1,mu2)(deltaBCIDistMetric(x1,x2,mu1,mu2));
else
	distancemetric = @(x,mu)((sum((x-repmat(mu,1,size(x,2))).^2,1)));
end
modelparams.distancemetric = distancemetric;
modelparams.counts = countmat;
tic
numreps = 50;
modelparams.numreps = numreps;
repeatfaflag = true;
while repeatfaflag
    LL = 0;
    allpershared= zeros(size(modelparams.counts,1),numreps);
    for n = 1:numreps
        display(['FA rep ',num2str(n)])
        [allestParams(n),tempLL] = fastfa(modelparams.counts ,zdim,'seed',n,'tol',1e-8);
        LL(n) = tempLL(end);
        allpershared(:,n) = diag(allestParams(n).L*allestParams(n).L')./(diag(allestParams(n).L*allestParams(n).L') + allestParams(n).Ph);
    end
    [~,bestLL]=max(LL);
    estParams = allestParams(bestLL);
    pershared = diag(estParams.L*estParams.L')./(diag(estParams.L*estParams.L') + estParams.Ph);
    repeatfaflag = false; % hopefully fewer FA dims fixes the high perc shared issue
%     if max(pershared)<maxpershared
%         repeatfaflag = 0;
%     else
%         modelparams.goodneuron(modelparams.goodneuron==1)= pershared<maxpershared;
%         modelparams.counts=modelparams.counts(pershared<maxpershared,:);
%         for n = 1:length(trainingdat)
%             trainingdat(n).counts = trainingdat(n).counts(pershared<maxpershared,:);
%         end
%     end
    figure
    errorbar(mean(allpershared,2),std(allpershared,[],2));
    hold on
    scatter(1:length(pershared),pershared,'MarkerEdgeColor','r')
    ylim([0 1])
end
goodneurons= modelparams.goodneuron;
toc
%% orthonormalize L and save evalues %%
L = estParams.L;
shared = L*L';
[evec,eval]=eig(shared);
deval = diag(eval);
tevec = evec*eval.^0.5;% orthogonalize (not orthonormalize) so that tevec*tevec'=L*L';
[sdeval,I]=sort(deval,'descend');
sevec = tevec(:,I);
estParams.L = sevec(:,1:zdim);
estParams.evals = sdeval;

L = estParams.L;
Ph = estParams.Ph;
sharedVar = diag(L*L');
percShared = sharedVar./(sharedVar+Ph) .* 100;
figure;
stem(percShared);
xlabel('channel #'); ylabel('% shared variance');
modelparams.estParams = estParams;
modelparams.zdim = zdim;
modelparams.trialnum = trialnum;
modelparamsold = modelparams;
%% create model
if ~isempty(dimsbci1)
    modelparams1 = computeDistanceModel(modelparamsold,trainingdat,dimsbci1,dimsbci2,deltabciflag);
    saveName = sprintf(['%s%s_lowdHmm'],savepath,filename);
    %modelparamsold = modelparams;
    modelparams1.bciind = 1;
    modelparams1.dimsbci = dimsbci1;
    modelparams = modelparams1;
else
    modelparams1 = [];
end

if ~isempty(dimsbci2)
    modelparams2 = computeDistanceModel(modelparamsold,trainingdat, dimsbci2,dimsbci1,deltabciflag);
    saveName = sprintf(['%s%s_lowdHmm'],savepath,filename);
    modelparams2.bciind = 2;
    modelparams2.dimsbci = dimsbci2;
    modelparams(2) = modelparams2; % if there are two bcis, this ensure that modelparams2 is second element
    
else
    modelparams2 = [];
end
%modelparams(1) = modelparams1;% this line works whether or not modelparams2 is added
save([saveName,num2str(trialnum),'.mat'],'modelparams');
end