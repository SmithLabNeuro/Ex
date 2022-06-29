function [trialstarts,trialends,trialstartinds,trialendinds]=detectMissingStartEndCode(trialstarts,trialends)
trialstarts = trialstarts(trialstarts<trialends(end));
mindistpair = zeros(length(trialstarts),1);minpairend = zeros(length(trialstarts),1);
for n = 1:length(trialstarts)
    codediffs = trialends-trialstarts(n);
    temp = find(codediffs>0,1);
    minpairend(n) = find(codediffs>0,1);
    mindistpair(n) = codediffs(minpairend(n));
end
a = unique(minpairend);
[minn,bin] = histc(minpairend,a);
multiple = find(minn>1);
goodtrials = ones(length(trialstarts),1);
for m = 1:length(multiple)
    b = find(ismember(bin, multiple(m))); 
    [~,c] = min(mindistpair(b));
    goodtrials(b((1:length(b))~=c))=0;
end

trialstarts = trialstarts(goodtrials==1);
trialstartinds = find(goodtrials);
trialends = trialends(a);
trialendinds = a;
end