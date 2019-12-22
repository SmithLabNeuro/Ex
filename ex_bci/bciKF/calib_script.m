function modelparams = calib_script(filename,filepath,nTrainTrials,nDelayBins,zDim,verbose,sortFlag,gamma)

    if ~exist('verbose','var'); verbose = true; end
    if ~exist('sortFlag','var'); sortFlag = 0; end
    if ~exist('gamma','var'); gamma = 0.2; end

    modelparams.nDelayBins = nDelayBins;
    modelparams.zDim = zDim;
    
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
    
    %% load data
    addpath(genpath('nev2bcistruct'));
    addpath(genpath('../structBuilders/'));
    addpath('../../../../smithlab/matlab/neural_net_sort_stuff');
    set(groot,'defaultAxesColorOrder',distinguishable_colors(16))
    
    %filepath = '../';
    if ~exist('filepath','var')||isempty(filepath)
        filepath = '../../../../../../../Volumes/DATA/wakko/';
    end
    if sortFlag == 1
        dat = nev2sortedStruct([filepath,filename,'.nev'],gamma,0);
    else
        nev = readNEV([filepath,filename,'.nev']);
        dat = nev2bcistruct(nev,'nevreadflag',1);
    end
    dat = dat([dat.result]==150);
    trainingdat = dat(1:nTrainTrials);
    
    bin = 0.05;
    fixateCode = 140;
    targOnCode = 70;
    targOffCode = 100;
    fixOffCode = 3;
    minchannelnum = 0;
    
    %% get angle shown on each trial
    for ii = 1:length(trainingdat)
        trainingdat(ii).angle = trainingdat(ii).params.trial.angle;
    end
    
    %% find good neurons (mean, fano, no cross-talk)
    if verbose
        fprintf('\tRemoving neurons with low FRs or crosstalk...\n');
    end
    
    addpath('goodchans');
    
    [~,tmp] = prepCalibCounts(trainingdat,fixateCode,fixOffCode,0,0.001);
    good_chans = get_good_channels(tmp);
    n_chans = length(good_chans);
    n_arr1 = sum(good_chans<=96);
    n_arr2 = sum(good_chans>96);
    modelparams.goodneuron = false(size(tmp(1).counts,1),1);
    modelparams.goodneuron(good_chans) = true;
    
	if verbose
        fprintf('\t\t%d channels remaining (%d in array 1, %d in array 2)\n',n_chans,n_arr1,n_arr2);
    end
    
    %% calibrate Kalman filter decoder (FA+KF)
    addpath('factor_analysis');
    
    % train fa model
    [~,fa_dat] = prepCalibCounts(trainingdat,targOnCode,fixOffCode,minchannelnum,bin);
    minBins = min([fa_dat.nBins]);
    x = nan(sum(modelparams.goodneuron),minBins,length(fa_dat));
    for i_trial = 1:length(fa_dat)
        x(:,:,i_trial) = fa_dat(i_trial).counts(modelparams.goodneuron,1:minBins);
    end
    
    x = reshape(x,size(x,1),[]);
    modelparams.fa_params = fastfa(x,zDim);
    
    % train kf model
    [~,kf_dat] = prepCalibCounts(trainingdat,targOffCode,fixOffCode,minchannelnum,bin);
    modelparams = train_lowd_kf(kf_dat,modelparams,'nDelayBins',nDelayBins);
    
    %% save model parameter file
    saveName = sprintf('X:satchel/%s_lowdKF.mat',filename(1:8));
    save(saveName,'modelparams');

end

