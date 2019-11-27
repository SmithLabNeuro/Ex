%% constants
% trial codes
correct = 150; bcicorrect = 161; saccade = 141;
filename = 'bcistruct_hmm';
numlatents = 10;
rewardCriteria = 0.25;
numintarg = 8;
%% load data
load(filename);
calibdata = bcistruct.calibrationdata;
bcidata = bcistruct.bcidata;
%% combine files 
for n = 1:length(calibdata); calibdata(n).calibtrial = 1;end % so we can toggle how train FA
saccadecodes = zeros(1,length(bcidata));
bcicorrectcodes = zeros(1,length(bcidata));
for n = 1:length(bcidata)
    bcidata(n).calibtrial = 0;
    saccadecodes(n) = double(~isempty(find(bcidata(n).trialcodes(:,2)==141,1)));
    bcicorrectcodes(n) = double(~isempty(find(bcidata(n).trialcodes(:,2)==161,1)));
end
% use correct trials only
correctcalib = calibdata([calibdata.result]==150);
% use bci_correct, bci missed w/ correct saccade, bci missed w/ incorrect saccade, and correct memory guided saccade
correctbci = bcidata(([bcidata.bcitrial]==0&[bcidata.result]==150)|...
                     ([bcidata.bcitrial]==1&(bcicorrectcodes==1|saccadecodes==1)));
bcifields=fieldnames(correctbci);calibfields=fieldnames(correctcalib);
for n = 1:length(bcifields);if ~isfield(correctcalib,bcifields{n}); correctcalib(1).(bcifields{n})=[];end;end
for n = 1:length(calibfields);if ~isfield(correctbci,calibfields{n}); correctbci(1).(calibfields{n})=[];end;end
dat = [correctcalib correctbci];
%% bin spikes in 50 ms bins
[~, ~, ~,dat,~,~] = prepCalibCounts(dat, 100, 3,0,0.05,10000);
%% compute psth for each condition
for n = 1:length(dat); dat(n).angle = dat(n).params.trial.angle;end
angle = [dat.angle];unangle = unique(angle);
for n = 1:length(unangle);psthstruct(n).psth = psth2017(dat(angle==unangle(n)));end
%% subtract psth from each trials spike counts (50 ms bins)
for n = 1:length(dat);dat(n).psthcounts=dat(n).counts-psthstruct(find(angle(n)==unangle,1)).psth(:,1:size(dat(n).counts,2));end
%% form spike count matrix and time matrix. maybe average across time within a trial
% average across bins in each trial (problem is unequal number of bins)
facounts = zeros(size(dat(1).psthcounts,1),length(dat));
for n = 1:length(dat); facounts(:,n) = mean(dat(n).psthcounts,2);end
pcacounts = zeros(size(dat(1).psthcounts,1),length(dat));
for n = 1:length(dat); pcacounts(:,n) = mean(dat(n).counts,2);end
for n = 1:length(dat); dat(n).recordingtime=dat(n).time(2)+dat(n).nevinfo.nevclockstart/30000;end
% throw everything into FA (problem is many zeros)
%% compute fa parameters
estparams = fastfa(facounts,numlatents);
%% extract top factor
L = estparams.L; shared = L*L'; [faevec,faeval] = eig(shared);
[sorteval,I] = sort(diag(faeval),'descend');fasortevenc=faevec(:,I);
topfactor = fasortevenc(:,1);
%% plot eigenspectrum
figure; plot(sorteval./sum(sorteval));xlabel('Eigen Index');ylabel('% Shared Var Expl');
figure; plot(sorteval./(sum(sorteval)+sum(estparams.Ph)));xlabel('Eigen Index');ylabel('% Total Var Expl');
%% project the data onto this factor
projection = topfactor'*facounts;
%% plot data
window = 60*10; stepsizemin = 60;
tic;
[smoothedbins,~,endsteps] = smoothSparse([dat.recordingtime],projection,0,stepsizemin, window,0,max([dat.recordingtime]));
toc;
figure;plot(endsteps/60,smoothedbins)

%% plot pca
pcatrainindices = [dat.calibtrial]==1|[dat.calibtrial]==0;
[coeff, score,latent] = pca(pcacounts(:,pcatrainindices)');
[smoothedbins,~,endsteps] = smoothSparse([dat(pcatrainindices).recordingtime],score(:,1)',0,stepsizemin, window,0,max([dat(pcatrainindices).recordingtime]));
figure;plot(endsteps/60,smoothedbins);
figure;plot(latent/sum(latent));
pcaparams.meanpcacounts =mean(pcacounts(:,pcatrainindices),2);
pcaparams.coeff = coeff;
%% plot population rates over time w/ and w/o pca 1
newscore = (coeff'*(pcacounts-repmat(pcaparams.meanpcacounts,1,size(pcacounts,2))))';
test = mean(coeff(:,2:end)*newscore(:,2:end)'+repmat(pcaparams.meanpcacounts,1,size(pcacounts,2)),1);
[smoothedbinstest,~,endsteps] = smoothSparse([dat.recordingtime],mean(test,1),0,stepsizemin, window,0,max([dat.recordingtime]));
[smoothedbinsall,~,~] = smoothSparse([dat.recordingtime],mean(pcacounts,1),0,stepsizemin, window,0,max([dat.recordingtime]));
figure;plot(endsteps,smoothedbinstest); hold on;plot(endsteps,smoothedbinsall);ylabel('Population FR');xlabel('Time (min)');
[smoothedbins,~,endsteps] = smoothSparse([dat.recordingtime],newscore(:,3)',0,stepsizemin, window,0,max([dat.recordingtime]));
figure;plot(endsteps/60,smoothedbins);
%% recompute bci hmm values
modelparams.allangles = unangle;
modelparams.goodneuron = find(ones(96,1));
datnodrift = dat;
for n = 1:length(datnodrift)
    lowdproj = (pcaparams.coeff(:,2:end)'*(dat(n).counts-repmat(pcaparams.meanpcacounts,1,size(dat(n).counts,2))))';
    datnodrift(n).counts =coeff(:,2:end)*lowdproj'+repmat(pcaparams.meanpcacounts,1,size(dat(n).counts,2));
end
modelparams = train_lowdHmm(dat([dat.calibtrial]==1),modelparams,'nDelayBins',6,'zDim',3,'hmmCovType','full');
datmapped = mapHMMBCI(dat([dat.calibtrial]==0),1:96,modelparams);
modelparamsnodrift = train_lowdHmm(datnodrift([dat.calibtrial]==1),modelparams,'nDelayBins',6,'zDim',3,'hmmCovType','full');
datnodriftmapped = mapHMMBCI(datnodrift([dat.calibtrial]==0),1:96,modelparamsnodrift);
bcidatmapped = datmapped([datmapped.bcitrial]==1);
bcidatnodriftmapped = datnodriftmapped([datnodriftmapped.bcitrial]==1);
for n = 1:length(bcidatmapped)
    rewardrange = bcidatmapped(n).bciStats.annulusProp<rewardCriteria;
    filterbank = ones(1,numintarg); bcidatmapped(n).bciresults = ~isempty(find(filter(filterbank,1,rewardrange)>=(numintarg),1));
end
for n = 1:length(bcidatnodriftmapped)
    rewardrange = bcidatnodriftmapped(n).bciStats.annulusProp<rewardCriteria;
    filterbank = ones(1,numintarg); bcidatnodriftmapped(n).bciresults = ~isempty(find(filter(filterbank,1,rewardrange)>=(numintarg),1));
end
plotSmoothBCICorrect(bcidatmapped)
plotSmoothBCICorrect(bcidatnodriftmapped)
