function modelparams = modelPNB(trainingdat,modelparams,traininds)
if ~exist('traininds','var'); traininds = 0; end
goodneuron = modelparams.goodneuron;
labels = [trainingdat.angle];
unlabels = trainingdat(1).params.block.delayAngle;
cellcounts = cell(length(unlabels),1);
for labelind = 1:length(unlabels)
    tempdat = trainingdat(labels==unlabels(labelind));
    tempcounts = [];
    for datind = 1:length(tempdat)
        if traininds ~= 0
            tempcounts = [tempcounts tempdat(datind).counts(:,traininds)];
        else
            tempcounts = [tempcounts tempdat(datind).counts];
        end
    end
    cellcounts{labelind} = tempcounts;
    % remove neurons with no spikes
    meancounts = find(mean(tempcounts,2)==0);
    if ~isempty(meancounts)
        goodneuron(meancounts) = 0;
    end
end

for labelind = 1:length(unlabels)
    modelparams.mean(:,labelind) = mean(cellcounts{labelind}(goodneuron==1,:),2);
end
modelparams.goodneuron = goodneuron;
end