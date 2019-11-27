function modelparams = computeDistanceModel(modelparams,trainingdat,dimsbci,dimsbci2,deltabciflag)


%% smooth distances MODEL SPECIFIC CODE
estParams = modelparams.estParams;
waitbins = modelparams.waitbins;
numbins = modelparams.numbins;

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
    if deltabciflag == 1
        smoothdist = modelparams.distancemetric(thissmoothlatent(dimsbci,:),thissmoothlatent(dimsbci2,:),modelparams.meanlatent(dimsbci),modelparams.meanlatent(dimsbci2));        
    else
        smoothdist = modelparams.distancemetric(thissmoothlatent(dimsbci,:),modelparams.meanlatent(dimsbci));
    end
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
difflist = correctlist-modelparams.targetCorrect;
[~,threshind]=min(abs(difflist));
goodthresh = threshold(threshind);
modelparams.threshold = goodthresh;
modelparams.allthresh = threshold;
modelparams.percentilevalues = values;
modelparams.percentileatthresh = values(threshind);

fprintf('Calib Percent Correct is %d \n',correctlist(threshind));

outdat = trainingdat;


%% determine offline correct (this is a sanity check and is not used in modelparams).
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

end