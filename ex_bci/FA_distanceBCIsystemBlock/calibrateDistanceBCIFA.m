function [modelparams, outdat] = calibrateDistanceBCIFA(filename,filepath,savepath,gamma,targetCorrect,multflag,offsetval,offlinemodeflag,trialnum,dat,modelparamsold,calibbin,numbins,waitbins,zdim,newcalcflag)
% need filepath
% gamma
% offlinemodeflag

% [modelparams, outdat] = calibrateDistanceBCIFA([],'D:/pepelepew/','X:/pepelepew/',0.2,0.5,0,0,0,0,[],[],[],[],[],5,removechans)
if ~exist('trialnum','var')|| isempty(trialnum); trialnum= 0;end
if ~exist('calibbin','var')|| isempty(calibbin);calibbin = 0.05;end
if ~exist('numbins','var')|| isempty(numbins);numbins = 8;end
if ~exist('waitbins','var')|| isempty(waitbins);waitbins = 8;end
if ~exist('zdim','var')|| isempty(zdim);zdim = 4;end
if ~exist('newcalcflag','var')|| isempty(newcalcflag);newcalcflag = 1;end
%filename = [];
if offlinemodeflag == 0
%% find file
if trialnum==0
if isempty(filename)
    thismonkey = filepath(4:5); %folder has format D:monkeyname
    thismonkey(1) = upper(thismonkey(1));
    thisyear = num2str(year(date));
    thisyear = thisyear(3:4);
    thismonth = num2str(month(date),'%02d');
    thisday = num2str(day(date),'%02d');
    thefile = dir([filepath,thismonkey,thisyear,thismonth,thisday,'*.nev']);
    
    if length(thefile)>1
        for n = 1:length(thefile)
            [~,thisfilename,~] = fileparts(thefile(n).name);
            userChoice = nan;
            while ~(strcmpi(userChoice,'y') || strcmpi(userChoice,'n'))
                userChoice = input([thisfilename,'.nev\nThis file? [y/n] ', ],'s');
            end
            
            switch lower(userChoice)
                case 'y'
                    filename = thisfilename;
                    break
                case 'n'
                otherwise
                    error('This shouldnt ever throw an error...something weird is going on');
            end
        end
        if isempty(filename)
            error('Please rerun code with desired filename as first input (do not include .nev extension)')
        end
    else
        
        [~,filename,~] = fileparts(thefile.name);
        fprintf(['Using the following file: \n',filename,'.nev\n'])  
        fprintf('If you want a different file, enter the filename\n as the first input to this function (do not include .nev extention)\n')
    end
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

countbinsize = calibbin;
modelparams.targetCorrect = targetCorrect;
modelparams.countbinsize = countbinsize;
modelparams.numbins = numbins;
modelparams.waitbins = waitbins;
waitbins = waitbins - 1; % 2018/11/28, RW do this b/c there is a communication delay between control and bci computer so the bci task only waits for waitbins-1 bins
modelparams.highscale = 5;
if trialnum == 0||newcalcflag ==0
    results = [dat.result];
    trainingdat = dat(results==150);
    onlythesetrials = 1:length(trainingdat);
    modelparams.onlythesetrials = onlythesetrials;
    trainingdat = trainingdat(onlythesetrials);
    [countmatpre, ~, ~,trainingdat,~,~] = prepCalibCounts(trainingdat,140, 3,0,countbinsize,1000);
    
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
    [~,~,~,tmpDat,~,~] = prepCalibCounts(trainingdat,140,3,0,0.001,[]);
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
    threshVals = meanVals+10*stdVals;
    
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
    countmat = [];
    modelparams.goodneuron = modelparamsold.goodneuron;
    modelparams.gamma = modelparamsold.gamma;
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
distancemetric = @(x,mu)((sum((x-repmat(mu,1,size(x,2))).^2,1)));
modelparams.distancemetric = distancemetric;
modelparams.counts = countmat;
maxpershared = 0.7;
tic
numreps = 50;
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
%% smooth distances MODEL SPECIFIC CODE
allsmoothdist = [];
allsmoothlatents = [];
for n = 1:length(trainingdat)
    thiscountunsmooth = trainingdat(n).counts;%(:,waitbins:end);
    z = fastfa_estep(thiscountunsmooth,estParams);
    trainingdat(n).thislatent = z.mean;
    smoothlatent = zeros(size(trainingdat(n).thislatent));
        for scoreind = 1:size(trainingdat(n).thislatent,2)
            if scoreind ==1
                smoothlatent(:,scoreind) = trainingdat(n).thislatent(:,scoreind);
            else
                smoothlatent(:,scoreind) = modelparams.alpha*trainingdat(n).thislatent(:,scoreind) +(1-modelparams.alpha)*smoothlatent(:,scoreind-1);
            end
        end
    trainingdat(n).smoothlatent =smoothlatent;
    allsmoothlatents = [allsmoothlatents smoothlatent];
end

modelparams.meanlatent = mean(allsmoothlatents,2);
for n = 1:length(trainingdat)
    thissmoothlatent = trainingdat(n).smoothlatent;
    smoothdist = distancemetric(thissmoothlatent,modelparams.meanlatent);
    trainingdat(n).smoothdist = smoothdist;
    allsmoothdist = [allsmoothdist smoothdist];
end


%% find threshold
stepsize = 0.1;
values = stepsize:stepsize:100;
threshold = zeros(length(values),1);
correctlist = zeros(length(values),1);
for  m = 1:length(threshold)
threshold(m) = prctile(allsmoothdist,values(m));
correct = zeros(length(trainingdat),1);
for n = 1:length(trainingdat)
    smoothdist = trainingdat(n).smoothdist(waitbins:end);
    correct(n) = ~isempty(find(filter(ones(1,numbins),1,(smoothdist)<threshold(m))==numbins,1));
end

correctlist(m) = (sum(correct)/length(correct));
end
difflist = correctlist-targetCorrect;
[~,threshind]=min(abs(difflist));
goodthresh = threshold(threshind);
modelparams.threshold = goodthresh;
modelparams.allthresh = threshold;
modelparams.percentilevalues = values;
modelparams.percentileatthresh = values(threshind);
modelparams.multgain = offsetval;
modelparams.muold = modelparams.mu;
if multflag == 1
modelparams.mu = modelparams.muold*modelparams.multgain;
else
    modelparams.mu = max(0,modelparams.muold+modelparams.multgain);
end
fprintf('Calib Percent Correct is %d \n',correctlist(threshind));
saveName = sprintf(['%s%s_lowdHmm'],savepath,filename);
save([saveName,num2str(trialnum),'.mat'],'modelparams');
outdat = trainingdat;


%% determine offline correct
for n = 1:length(trainingdat)
    thissmooth = trainingdat(n).smoothdist(waitbins:end);
    thispercent = zeros(size(thissmooth));
    for m = 1:length(thispercent)
    thispercent(m) = max(1,sum(modelparams.allthresh<thissmooth(m)));
    end
    valtosend = modelparams.percentilevalues(thispercent)/modelparams.percentileatthresh;
    valtosend2 = zeros(size(valtosend));
    valtosend3 = zeros(size(valtosend));
    for m = 1:length(valtosend)
    if valtosend(m) < 1
        valtosend2(m) = (valtosend(m)-1)*modelparams.percentileatthresh/(modelparams.percentileatthresh-10)+1;
    else
        valtosend2(m) = (valtosend(m)-1)*(modelparams.highscale-1)*modelparams.percentileatthresh/(90-modelparams.percentileatthresh)+1;
    end
        valtosend3(m) = max(0,min(1,0.2*valtosend2(m)));
    end
    
    endtrialbin = find(filter(ones(modelparams.numbins,1),1,valtosend3<0.2)==8,1,'first');
    if ~isempty(endtrialbin)
        trainingdat(n).offlinecorrect = 1;
    else
        trainingdat(n).offlinecorrect = 0;
    end
    trainingdat(n).valtosend = valtosend3;
end

% figure
% subplot(3,1,1)
% plot(trainingdat(10).valtosend)
% subplot(3,1,2)
% plot(trainingdat(20).valtosend)
% subplot(3,1,3)
% plot(trainingdat(30).valtosend)

end