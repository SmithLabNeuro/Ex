function [neuronstokeep, dirsel, correlations] = dirselNeuronToKeep(dat,varargin)
% corr(even_means odd_means) >0.25, separate out conditions first, then
% assign trials to even or odd groups to ensure equal numbers in each group
% 300 ms from targ off to min delay end.
% load Pe180712_s374a_dirmem_bci_discrete_0001_sort0_2

%% make sure important stuff in path
addpath(genpath('../structbuilders'))

%% optional parameters
channels = [1:96 257:352];
startcode = 100;
endcode = 3;
startoffset = 0.3;
endoffset = 0;
mincorr = 0.25;
sortflag = 1;
minrate = 2;
prctThresh = 90;
dirselflag = false;
goodratesflag = true;
goodcorrelationsflag = false;
goodDepthFlag = true;
assignopts(who,varargin);

%% make sure only looking at correct trials
results = [dat.result];
gooddat = dat(results==150);

%% get conditions
angles = driftchoiceextractparam(gooddat,'angle=');
unangles = unique(angles);

%% get min end time
trialstarts = zeros(length(gooddat),1);
endtimes = zeros(length(gooddat),1);
for n = 1:length(gooddat)
    trialcodes = gooddat(n).trialcodes;
    trialstarts(n) = trialcodes(find(trialcodes(:,2)==startcode,1),3)+startoffset; % align to fix time
    endtimes(n) = trialcodes(trialcodes(:,2)==endcode,3)-trialstarts(n);
end
endtime = min(endtimes) + endoffset;
%% compute counts
allcounts = zeros(length(channels),length(gooddat));
for n = 1:length(gooddat)
    spikes = unpackSpikes(gooddat(n).spikeinfo, gooddat(n).spiketimesdiff,gooddat(n).firstspike,30000);
    thisspikes = spikes(spikes(:,3)>trialstarts(n) & spikes(:,3)<(trialstarts(n)+endtime),:);
    tempcounts = histc(thisspikes(thisspikes(:,2)==sortflag,1),channels);
    allcounts(:,n) = tempcounts;
end

%% compute means
odd = zeros(length(channels), length(unangles));
even = zeros(length(channels), length(unangles));
alltune = zeros(length(channels), length(unangles));
for n = 1:length(unangles)
    angleinds = find(angles==unangles(n));
    temp = 1:length(angleinds);
    oddinds = angleinds(mod(temp,2)==1);
    eveninds = angleinds(mod(temp,2)==0);
    if length(oddinds)>length(eveninds)
        oddinds = oddinds(1:end-1);
    end
    odd(:,n) = mean(allcounts(:,oddinds),2);
    even(:,n) = mean(allcounts(:,eveninds),2);
    alltune(:,n) = mean(allcounts(:,angleinds),2);
end

%% compute correlations
correlations = zeros(length(channels),1);
dirsel = zeros(length(channels),1);
for n = 1:length(correlations)
    correlations(n) = corr(odd(n,:)',even(n,:)');
    dirsel(n) = orivecfit(unangles, alltune(n,:), 0);
end

%% compute permutation test for dirsel and modulation depth
numiters = 1000;
randdirsel = zeros(size(alltune,1),numiters);
alltunerand = zeros(size(alltune));
for n = 1:numiters   
    randangles = angles(randperm(length(angles)));
    for m = 1:length(unangles)
        alltunerand(:,m) = mean(allcounts(:,randangles==unangles(m)),2);
    end
    for m = 1:size(alltunerand,1)
        randdirsel(m,n) = orivecfit(unangles, alltunerand(m,:), 0);
    end
end
if dirselflag
    sigdirsel = dirsel>prctile(randdirsel,prctThresh,2);
else
    sigdirsel = ones(size(dirsel))==1;
end
if goodratesflag
    goodrates = (mean(allcounts,2)/(endtime-startoffset))>minrate;
else
    goodrates = ones(size(mean(allcounts,2)))==1;
end
if goodcorrelationsflag
    goodcorrelations = correlations>mincorr;
else
    goodrates = ones(size(mean(allcounts,2)))==1;
end
neuronstokeep = (sigdirsel&goodrates&goodcorrelations);
end