function percentcorrect = evalPNBmodel(dat,modelparams)
means = modelparams.mean;
goodneuron = modelparams.goodneuron;
angles = dat(1).params.block.delayAngle;
mintrials = 10;
% determine max number of elements 
maxel = 0;
for n = 1:length(dat)
    curel = size(dat(n).counts,2);
    if curel>maxel
        maxel = curel;
    end
end
predictedlabels = nan(length(dat),maxel);
actuallabels = nan(length(dat),maxel);
correctmat = nan(length(dat),maxel);
for n = 1:length(dat)
    counts = dat(n).counts(goodneuron,:);
    thisind = find(dat(n).angle==angles);
    for m = 1:size(counts,2)
        [~,templabel] = max(computePostNBP(counts(:,m),means));
        predictedlabels(n,m) = templabel;
        actuallabels(n,m) = thisind;
        correctmat(n,m) = templabel == thisind;
    end
end
numnonnan = sum(~isnan(correctmat),1);
percentcorrect = nansum(correctmat(:,numnonnan>mintrials),1)./numnonnan(numnonnan>mintrials);

end