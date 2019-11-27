function modelparams = calibBCI_online_discrete(filename,filepath,cross_channel_flag,ind1,ind2,numtrials)
%% load data
% if offline mode data, separate into train and test data.
addpath(genpath('nev2akashstruct'));
addpath ../
%filepath = '../';
if ~exist('filepath','var')||isempty(filepath)
    filepath = '../../../../../../../Volumes/DATA/wakko/';
end
if ~exist('ind1','var') || isempty(ind1); ind1 = 1; end
if ~exist('ind2','var') || isempty(ind2); ind2 = 0; end

%21 119 40 (green color in bci)
%filename = 'Wa171102_s169a_dirmem_withhelp_0001';
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
dat = nev2akashstruct([],nev,1,1,0);
results = [dat.result];
correctdat = dat(results==150);

angletemp = driftchoiceextractparam(correctdat,'angle=');
modelparams.allangles = dat(1).params.block.delayAngle;
% correctdat = correctdat(angletemp==45 | angletemp==135 |angletemp==225 | angletemp==315);
% correctdat(1).params.block.delayAngle = [45 135 225 315];
othername = '50msalltrialonlinetest';
bin = 0.4;
targetoff = 100;
mintimetotarget = 1;
fixoff = 3;
minchannelnum = 0;
% convert to x,y?
if exist('numtrials','var'); 
    trainingdat = correctdat(1:numtrials);
    if numtrials < length(correctdat)
        testdat = correctdat(numtrials+1:end);
    end
else
    trainingdat = correctdat;
end


%% define data (Make this into a function so that you can reuse for bci dat)
% NOTE THIS LIMITS SPIKE BINS TO BE 5 bins long per trial
[countmat, ~, ~,trainingdat,~,~] = prepCalibCounts(trainingdat, 70, fixoff,minchannelnum,bin,100);


%% calibrate decoder (poisson)
modelparams.goodneuron =mean(countmat,2)>(1*bin) & (var(countmat,[],2)./mean(countmat,2))<5;% 
modelparams.goodneuron(97:end) = 0;
%% calibrate kalman filter
%% get projections 
%[~, ~,trainingdat] = positionStates(trainingdat,mintimetotarget/bin,1,1);
if ~exist('numtrials','var'); 
    modelparams = modelPNB(trainingdat,modelparams,ind1:ind2);
    percentcorrect = evalPNBmodel(trainingdat,modelparams);
    mean(percentcorrect)
else
    modelparams = modelPNB(trainingdat,modelparams,1:numtrials);
    percentcorrect = evalPNBmodel(trainingdat(1:numtrials),modelparams);
    mean(percentcorrect)
end

save([filename,'modelparamstestgauss',othername],'modelparams');

% need to add decoding accuracy output
%traindatnew = dat2Traj(trainingdat,modelparams);
if exist('numtrials','var') && numtrials < length(correctdat)
    [~, ~, ~,testdat,~,~] = prepCalibCounts(testdat, targetoff, fixoff,minchannelnum,bin,100);
    percentcorrect_test = evalPNBmodel(testdat,modelparams)
    %testdat = dat2Traj(testdat,modelparams);
end

targofftime = zeros(length(trainingdat),1);
fixofftime = zeros(length(trainingdat),1);
for n = 1:length(trainingdat)
    targontime = trainingdat(n).events(trainingdat(n).events(:,2)==70,3);
    thistargofftime = trainingdat(n).events(trainingdat(n).events(:,2)==targetoff,3);
    thisfixofftime = trainingdat(n).events(trainingdat(n).events(:,2)==fixoff,3);
    targofftime(n) = thistargofftime(1) - targontime(1);
    fixofftime(n) = thisfixofftime(1) - targontime(1);
end

figure
scatter(bin/2+(0:bin:(length(percentcorrect)-1)*bin),percentcorrect,'FaceColor','k');
hold on
ylim([0 1])
line([mean(targofftime) mean(targofftime)], ylim,'Color','r');
line([mean(fixofftime) mean(fixofftime)], ylim, 'Color','b');

end