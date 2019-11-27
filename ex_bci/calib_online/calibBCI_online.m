function modelparams = calibBCI_online(filename,filepath,cross_channel_flag,ind1,ind2,numtrials)
%% load data
% if offline mode data, separate into train and test data.
addpath(genpath('nev2akashstruct'));
%filepath = '../';
if ~exist('filepath','var')||isempty(filepath)
    filepath = '../../../../../../../Volumes/DATA/wakko/';
end
if ~exist('ind1','var') || isempty(ind1); ind1 = 1; end
if ~exist('ind2','var') || isempty(ind2); ind2 = 0; end

%filename = 'Wa171102_s169a_dirmem_withhelp_0001';
nev = readNEV([filepath,filename,'.nev']);
if cross_channel_flag == 1
    digcodes = nev(nev(:,1)==0,:);
    spikes = nev(nev(:,1)~=0,:);
    timebin = 0.0005;
    numevents = 20;
    [~, percentexclude,excludecell,~] = excludeCrossChannelNoise(spikes,timebin,numevents,0);
    fprintf('Percent of spikes excluded: %g\n',percentexclude)
    spikes2 = spikes(excludecell==1,:);
    nev = [digcodes;spikes2];
    [~,I] = sort(nev(:,3),'ascend');
    nev = nev(I,:);
end
dat = nev2akashstruct([],nev,1,1,0);
results = [dat.result];
correctdat = dat(results==150);
distancetemp = driftchoiceextractparam(correctdat,'distance=');
correctdat = correctdat(distancetemp>80);
if ~exist('numtrials','var'); numtrials = length(correctdat);end
othername = '50msalltrialonlinetest';
bin = 0.05;
targetoff = 100;
mintimetotarget = 1;
fixoff = 3;
minchannelnum = 256;
% convert to x,y?
trainingdat = correctdat(1:numtrials);
if numtrials < length(correctdat)
    testdat = correctdat(numtrials+1:end);
end


%% define data (Make this into a function so that you can reuse for bci dat)
% NOTE THIS LIMITS SPIKE BINS TO BE 5 bins long per trial
[countmat, ~, ~,trainingdat,~,~] = prepCalibCounts(trainingdat, targetoff, fixoff,minchannelnum,bin,30);


%% calibrate decoder (poisson)
modelparams.goodneuron =mean(countmat,2)>(1*bin) & (var(countmat,[],2)./mean(countmat,2))<5;% 

%% calibrate kalman filter
%% get projections 
[~, ~,trainingdat] = positionStates(trainingdat,mintimetotarget/bin,1,1);
modelparams.KFparamsdec = ps9kftrain(trainingdat,modelparams,10,[],ind1,ind2);
save([filename,'modelparamstestgauss',othername],'modelparams');

traindatnew = dat2Traj(trainingdat,modelparams);
if numtrials < length(correctdat)
    [~, ~, ~,testdat,~,~] = prepCalibCounts(testdat, targetoff, fixoff,minchannelnum,bin,30);
    [~, ~,testdat] = positionStates(testdat,mintimetotarget/bin,1,1);
testdat = dat2Traj(testdat,modelparams);
end
end