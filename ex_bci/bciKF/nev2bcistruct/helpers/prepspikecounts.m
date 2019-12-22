function alltrials = prepspikecounts(dat,binsize,offset,channelinds,startalign, endalign,result)
% Parameters
if ~exist('binsize','var')
    binsize = 0.1;
end

keeptrials = zeros(length(dat),1);
if ~isempty(result)
for n = 1:length(dat)
    for m = 1:length(result)
        if ~isempty(find(dat(n).result ==result(m),1))
        keeptrials(n) = 1;
        end
    end
end
dat = dat(keeptrials==1);
end
% Constants
correct = 150;

keeptrials = ones(length(dat),1);

trialresults = [dat.result];

initcue = driftchoiceextractparam(dat(1),'initcue=');
miniblocksize = driftchoiceextractparam(dat(1),'miniblocksize=');
correctresults = cumsum(trialresults == correct);
correctresults(correctresults==0) = 1;
if initcue == 1
    cueinout = (1*(mod((floor((correctresults-1)/miniblocksize)+1),2)==1)+2*(mod((floor((correctresults-1)/miniblocksize)+1),2)==0))';
else
    cueinout = 3-(1*(mod((floor((correctresults-1)/miniblocksize)+1),2)==1)+2*(mod((floor((correctresults-1)/miniblocksize)+1),2)==0))';
end
correcttrials = dat(keeptrials==1);
cueinout = cueinout(keeptrials==1);

% get spike counts in 100 ms bins starting at stim onset
isvalid = driftchoiceextractparam(correcttrials,'isValid=');
oripicktemp = driftchoiceextractparam(correcttrials,'oriPick=');
oripick = convertOri(oripicktemp, isvalid, cueinout);

% block cueing
alltrials = driftchoicebinspikes(correcttrials,startalign,endalign,binsize,channelinds,offset);
%% Subtract PSTH

for n = 1:length(alltrials)  
    alltrials(n).cue = cueinout(n);
    alltrials(n).ori = oripick(n);
    alltrials(n).isvalid = isvalid(n);
    alltrials(n).correctind = correctresults(n);
end


end