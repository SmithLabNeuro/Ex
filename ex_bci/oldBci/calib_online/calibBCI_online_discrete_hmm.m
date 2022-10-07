function modelparams = calibBCI_online_discrete_hmm(filename,filepath,cross_channel_flag,nDelayBins,zDim,sortFlag,gamma,dirselflag)

if ~exist('sortFlag','var');sortFlag = 0;end
if ~exist('gamma','var');gamma = 0.2;end


%% find file
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
    mySavePath = ['C:/Users/smithlab/dropbox/smithlabdata/bcidailyplots/' subjectDate];
    if exist(mySavePath,'dir')
        warning([mySavePath ' folder already exists'])
    else
        mkdir(mySavePath)
    end
    
    %% load data
    % if offline mode data, separate into train and test data.
    addpath(genpath('../structBuilders/'));
    addpath('../bciHMM');
    addpath('../dailyPlots/dailyDiagnostics/bayesClassifier/');
    addpath('../');
    addpath('../../../../smithlab/matlab/neural_net_sort_stuff');
    addpath('../distinguishable_colors/');
    set(groot,'defaultAxesColorOrder',distinguishable_colors(16))
    
    %filepath = '../';
    if ~exist('filepath','var')||isempty(filepath)
        filepath = '../../../../../../../Volumes/DATA/wakko/';
    end
    if sortFlag == 1
        dat = nev2sortedStruct([filepath,filename,'.nev'],gamma,0);
    else
        nev = readNEV([filepath,filename,'.nev']);
        if cross_channel_flag == 1
            digcodes = nev(nev(:,1)==0,:);
            spikes = nev(nev(:,1)~=0,:);
            timebin = 0.0005;
            numevents = 10;
            [~, percentexclude,excludecell,~] = excludeCrossChannelNoise(spikes,timebin,numevents,0);
            fprintf('Percent of spikes excluded: %g\n',percentexclude)
            spikes2 = spikes(excludecell==1,:);
            nev = [digcodes;spikes2];
            [~,I] = sort(nev(:,3),'ascend');
            nev = nev(I,:);
        end
        dat = nev2bcistruct(nev,'nevreadflag',1);
    end
    trainingdat = dat([dat.result]==150);
    angletemp = driftchoiceextractparam(trainingdat,'angle=');
    modelparams.allangles = unique(angletemp);
    
    bin = 0.05;
    fixateCode = 140;
    targOffCode = 100;
    fixOffCode = 3;
    minchannelnum = 0;
    
    
    %% find good neurons
    % uniqueAngles = unique(angletemp);
    % [~, ~, ~,counts_2s,~,~] = prepCalibCounts(trainingdat,fixateCode,fixOffCode,minchannelnum,0.1,100);
    % allCounts = nan(size(counts_2s(1).counts,1),10,length(counts_2s));
    % for ii = 1:length(counts_2s)
    %     allCounts(:,:,ii) = counts_2s(ii).counts(:,1:10);
    % end
    % for i_ang = 1:length(uniqueAngles)
    %     currDat = allCounts(:,:,uniqueAngles(i_ang)==angletemp);
    %     fanoFact_cond(:,:,i_ang) = var(currDat,[],3)./mean(currDat,3);
    % end
    % meanFR = mean(allCounts,3)./0.1;
    % avgMeanFR = mean(meanFR,2);
    % goodMean = avgMeanFR>2 & avgMeanFR<200;
    % avgFano = mean(mean(fanoFact_cond,3),2);
    % goodFano = avgFano>.5 & avgFano<100;
    % modelparams.goodneuron = goodMean & goodFano;
    %
    % % 7/18/2018 RW added dirsel neuron selection criteria
    % if dirselflag == 1
    %     [neuronstokeep, dirsel, correlations] = neuronToKeep(trainingdat);
    %     modelparams.goodneuron = (neuronstokeep==1)&(modelparams.goodneuron==1);
    % end
    %% format training data
    [~, ~, ~,trainingdat,~,~] = prepCalibCounts(trainingdat, targOffCode, fixOffCode,minchannelnum,bin,100);
    
    %% find neurons that have good firing rates and significant modulation depths
    fprintf('Figuring out which neurons are good...\n');
    modelparams.goodneuron = neuronToKeep(trainingdat,'nDelayBins',nDelayBins);
    
    %% choose which targets to include
    % targAngs = [trainingdat.angle]';
    % keepAngles = unique(targAngs);
    % keepAngles = setdiff(keepAngles,135:22.5:225);
    % trainingdat = trainingdat(ismember(targAngs,keepAngles));
    
    %% decide which array(s) to use (bayes classifiers and confusion matrices)
    % train models
    [~,allPred.pred_both] = train_lowdBayes(trainingdat,modelparams.goodneuron,'whichArray','both','nDelayBins',nDelayBins,'zDim',zDim,'covType','sharedgauss');
    if sum(modelparams.goodneuron(97:end))>zDim
        [neuronsKept_left,allPred.pred_left] = train_lowdBayes(trainingdat,modelparams.goodneuron,'whichArray','left','nDelayBins',nDelayBins,'zDim',zDim,'covType','sharedgauss');
    else
        neuronsKept_left = modelparams.goodneuron;
        neuronsKept_left(97:end) = false;
        allPred.pred_left.y_true = allPred.pred_both.y_true;
        allPred.pred_left.y_cv = nan(length(allPred.pred_both.y_true),1);
        allPred.pred_left.y_training = nan(length(allPred.pred_both.y_true),1);
    end
    if sum(modelparams.goodneuron(1:96))>zDim
        [neuronsKept_right,allPred.pred_right] = train_lowdBayes(trainingdat,modelparams.goodneuron,'whichArray','right','nDelayBins',nDelayBins,'zDim',zDim,'covType','sharedgauss');
    else
        neuronsKept_right = modelparams.goodneuron;
        neuronsKept_right(1:96) = false;
        allPred.pred_right.y_true = allPred.pred_both.y_true;
        allPred.pred_right.y_cv = nan(length(allPred.pred_both.y_true),1);
        allPred.pred_right.y_training = nan(length(allPred.pred_both.y_true),1);
    end
    
    % plot confusion matrices for training predictions
    plotConfMat(allPred,'training');
    print([mySavePath '/arrayDecodingQuality_training.png'],'-dpng');
    
    % plot confusion matrices for held-out predictions
    plotConfMat(allPred,'crossval');
    print([mySavePath '/arrayDecodingQuality_crossval.png'],'-dpng');
    
    userChoice = nan;
    while ~(strcmpi(userChoice,'both') || strcmpi(userChoice,'left') || strcmpi(userChoice,'right'))
        userChoice = input('Which array should we use today (both,left,right)? ','s');
    end
    
    switch lower(userChoice)
        case 'both'
        case 'left'
            modelparams.goodneuron = neuronsKept_left;
        case 'right'
            modelparams.goodneuron = neuronsKept_right;
        otherwise
            error('Specify a valid array to use (both,left,right)');
    end
    
    %% calibrate hmm decoder (pLDA+HMM)
    modelparams = train_lowdHmm(trainingdat,modelparams,'nDelayBins',nDelayBins,'zDim',zDim,'hmmCovType','shareddiag');
    print([mySavePath '/calibration_LDA.png'],'-dpng');
    
    %% save model parameter file
    %saveName = sprintf('../../../../../../../Volumes/bciData/wakko/%s_lowdHmm.mat',filename);
    saveName = sprintf('X:pepelepew/%s_lowdHmm.mat',filename);
    save(saveName,'modelparams');
end