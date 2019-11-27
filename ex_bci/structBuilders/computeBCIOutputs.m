function dat = computeBCIOutputs(dat,modelparams,calcCorrectFlag,useinputlatent)
if ~exist('useinputlatent','var');useinputlatent=0;end
if ~exist('calcCorrectFlag','var')||isempty(calcCorrectFlag);calcCorrectFlag=1;end
allsmoothdist = [];
allsmoothlatents = [];
if isfield(modelparams,'dimsbci')
    dimsbci = modelparams.dimsbci;
else
    dimsbci = 1:length(modelparams.meanlatent);
end
for n = 1:length(dat)
    if useinputlatent == 0
        thiscountunsmooth = dat(n).counts;%(:,waitbins:end);
        z = fastfa_estep(thiscountunsmooth,modelparams.estParams);
        dat(n).thislatent = z.mean;
    end
    smoothlatent = zeros(size(dat(n).thislatent));
        for scoreind = 1:size(dat(n).thislatent,2)
            if scoreind ==1
                smoothlatent(:,scoreind) = dat(n).thislatent(:,scoreind);
            else
                smoothlatent(:,scoreind) = modelparams.alpha*dat(n).thislatent(:,scoreind) +(1-modelparams.alpha)*smoothlatent(:,scoreind-1);
            end
        end
    dat(n).smoothlatent =smoothlatent;
    dat(n).smoothdist = modelparams.distancemetric(smoothlatent(dimsbci,:),modelparams.meanlatent(dimsbci));
end
if calcCorrectFlag == 1
    for n = 1:length(dat)
        thissmooth = dat(n).smoothdist(1:end);
        thispercent = zeros(size(thissmooth));
        for m = 1:length(thispercent)
        thispercent(m) = max(1,sum(modelparams.allthresh<thissmooth(m)));
        end
        
        dat(n).percentile = modelparams.percentilevalues(thispercent);
        valtosend = modelparams.percentilevalues(thispercent)/modelparams.percentileatthresh;
        valtosend2 = zeros(size(valtosend));
        valtosend3 = zeros(size(valtosend));
        for m = 1:length(valtosend)
        if valtosend(m) < 1
            valtosend2(m) = (valtosend(m)-1)*modelparams.percentileatthresh/(modelparams.percentileatthresh-10)+1;
        else
            valtosend2(m) = (valtosend(m)-1)*(modelparams.highscale-1)*modelparams.percentileatthresh/(90-modelparams.percentileatthresh)+1;
        end
            valtosend3(m) = max(0,min(1,valtosend2(m)/modelparams.highscale));
        end
% Index valtosend3 starting at waitbins-1 to account for communication delay b/t control and bci. This should also be included in the calibration
% Index valtosend3 end-1 to account for a delay between the task ending andthe bci stops running
        endtrialbin = find(filter(ones(modelparams.numbins,1),1,valtosend3((modelparams.waitbins-1):(end-1))<=(1/modelparams.highscale))==8,1,'first');
        if ~isempty(endtrialbin)
            endtrialbin = endtrialbin + modelparams.waitbins-1;
        end
        if ~isempty(endtrialbin)
            dat(n).offlinecorrect = 1;
        else
            dat(n).offlinecorrect = 0;
        end
        dat(n).endtrialbin = endtrialbin;
        dat(n).valtosend = valtosend3;
    end
end
end