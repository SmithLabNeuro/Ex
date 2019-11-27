function [modelparams, outdat] = calibrateDistanceBCIgausecounts(filename,filepath,savepath,gamma,targetCorrect,multflag,offsetval,offlinemodeflag,trialnum,dat,modelparamsold,calibbin,numbins,waitbins)
% need filepath
% gamma
% offlinemodeflag
if ~exist('trialnum','var'); trialnum= 0;end
if ~exist('calibbin','var');calibbin = 0.05;end
if ~exist('numbins','var');numbins = 8;end
if ~exist('waitbins','var');waitbins = 2;end
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

modelparams.filtertimeconstanst = 4;
modelparams.alpha = 1-exp(-1/modelparams.filtertimeconstanst);
modelparams.targetCorrect = targetCorrect;
modelparams.countbinsize = countbinsize;
modelparams.numbins = numbins;
modelparams.waitbins = waitbins;
modelparams.highscale = 5;
if trialnum == 0
    results = [dat.result];
    trainingdat = dat(results==150);
    onlythesetrials = 1:60;%length(trainingdat);
    modelparams.onlythesetrials = onlythesetrials;
    trainingdat = trainingdat(onlythesetrials);
    [countmatpre, ~, ~,trainingdat,~,~] = prepCalibCounts(trainingdat,140, 3,0,countbinsize,1000);
    
    meanRate = mean(countmatpre,2).*(1/countbinsize);
    fanoFact = var(countmatpre,[],2)./mean(countmatpre,2);
    goodneuron = meanRate > 2 & fanoFact<3 & fanoFact>0.5;%meanRate > 2 & fanoFact<3 & fanoFact>0.5;
    goodneuron(97:end) = false;
    fprintf('   %d of %d are good neurons\n',sum(goodneuron),96);
    
    countmat = countmatpre(goodneuron,:);
    modelparams.goodneuron = goodneuron;
    modelparams.gamma = gamma;
    for n = 1:length(trainingdat)
        trainingdat(n).counts = trainingdat(n).counts(modelparams.goodneuron,:);
    end
    
else
    countmat = [];
    for n = 1:length(trainingdat)
        countmat = [countmat trainingdat(n).counts];
    end
    modelparams.goodneuron = modelparamsold.goodneuron;
    modelparams.gamma = modelparamsold.gamma;
    goodneuron = modelparams.goodneuron;%ones(size(countmat,1),1)==1;
end

%% define mu
modelparams.mu = mean(countmat,2);
distancemetric = @(x,mu)((sum((x-repmat(mu,1,size(x,2))).^2,1)));
modelparams.distancemetric = distancemetric;
modelparams.counts = countmat;
%% define kfparams %%%% NEED TO DECIDE HOW TO DEFINE PARAMETERS, 
KFparams.mutmin = 0;
KFparams.sigtmin = 0.1;
KFparams.A = 1;
KFparams.Q = 1;
KFparams.C(1) = 1;
KFparams.C(2) = 0;
KFparams.R = 2;
modelparams.KFparams = KFparams;


%% smooth distances
allsmoothdist = [];
for n = 1:length(trainingdat)
    thiscountunsmooth = trainingdat(n).counts;
    meanvar.mu = KFparams.mutmin;
    meanvar.sig = KFparams.sigtmin;
    thiscount = zeros(size(thiscountunsmooth));
    for scoreind = 1:size(thiscountunsmooth,2)
        for neurind = 1:size(thiscountunsmooth,1)
%         if scoreind ==1
%             smoothdist(scoreind) = distances(scoreind);
%         else
%             smoothdist(scoreind) = modelparams.alpha*distances(scoreind) +(1-modelparams.alpha)*smoothdist(scoreind-1);
%         end
        [meanvar] = smoothCounts(thiscountunsmooth(neurind,scoreind),KFparams,meanvar);
        thiscount(neurind,scoreind) = meanvar.mu;
        end
    end
    
    
    
    smoothdist = distancemetric(thiscount,modelparams.mu);
    
    % apply smoothing
%     meanvar.mu = KFparams.mutmin;
%     meanvar.sig = KFparams.sigtmin;
%     smoothdist = zeros(length(distances),1);
%     for scoreind = 1:length(distances)
%         if scoreind ==1
%             smoothdist(scoreind) = distances(scoreind);
%         else
%             smoothdist(scoreind) = modelparams.alpha*distances(scoreind) +(1-modelparams.alpha)*smoothdist(scoreind-1);
%         end
% %         [meanvar] = smoothCounts(distances(scoreind),KFparams,meanvar);
% %         smoothdist(scoreind) = meanvar.mu;
%     end
    trainingdat(n).smoothcounts = thiscount;
    trainingdat(n).smoothdist = smoothdist;
    allsmoothdist = [allsmoothdist smoothdist];
end

%% find threshold
stepsize = 1;
values = stepsize:stepsize:100;
threshold = zeros(length(values),1);
correctlist = zeros(length(values),1);
for  m = 1:length(threshold)
    m
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
save([saveName,num2str(trialnum),'_gain1_05.mat'],'modelparams');
outdat = trainingdat;


%% determine offline correct
for n = 1:length(trainingdat)
    thissmooth = trainingdat(n).smoothdist;
    thispercent = zeros(size(thissmooth));
    for m = 1:length(thispercent)
    thispercent(m) = max(1,sum(modelparams.allthresh<thissmooth(m)));
    end
    valtosend = thispercent/modelparams.percentileatthresh;
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
    
    endtrialbin = find(filter(ones(modelparams.numbins,1),1,valtosend2(modelparams.waitbins:end)<0.2)==8,1,'first');
    if ~isempty(endtrialbin)
        trainingdat(n).offlinecorrect = 1;
    else
        trainingdat(n).offlinecorrect = 0;
    end
    trainingdat(n).valtosend = valtosend3;
end

figure
subplot(3,1,1)
plot(trainingdat(10).valtosend)
subplot(3,1,2)
plot(trainingdat(20).valtosend)
% subplot(3,1,3)
% plot(trainingdat(30).valtosend)

end