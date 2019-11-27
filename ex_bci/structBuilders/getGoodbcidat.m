function bcidat = getGoodbcidat(bcidatpre,modelparams)
if ~exist('modelparams','var'); modelparams = [];end
results = [bcidatpre.result];
bcidat = bcidatpre(results==161 |results==162);
emptylog = zeros(length(bcidat),1);
for n = 1:length(bcidat)
    emptylog(n) = isempty(bcidat(n).trialstarttime);
end
bcidat = bcidat(emptylog==0);
bcidat = computeCountFromLog(bcidat);

if ~isempty(modelparams)
    for n = 1:length(bcidat)
        bcidat(n).counts = bcidat(n).logcount(modelparams.goodneuron,:);
    end
end
end