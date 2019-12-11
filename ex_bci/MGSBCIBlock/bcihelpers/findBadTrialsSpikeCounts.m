function [bcidat,trialstoremove]=findBadTrialsSpikeCounts(bcidat,threshold)
if ~exist('threshold','var');threshold = 6;end
sumcounts = [];
trialnums = [];
maxcount = [];
for n = 1:length(bcidat)
    sumcounts = [sumcounts, sum(bcidat(n).counts,1)];
    trialnums = [trialnums, n*ones(1,size(bcidat(n).counts,2))];
    maxcount = [maxcount, max(sum(bcidat(n).counts,1))];
end
meansumcounts = mean(sumcounts);
stdsumcounts = std(sumcounts);
minbin = min(sumcounts);
maxbin = max(sumcounts);
histbins=((minbin:maxbin)-meansumcounts)/stdsumcounts;
zscore = (sumcounts-meansumcounts)/stdsumcounts;
%numbins = length(unique(zscore));
figure
subplot(3,1,1)
[a,b]=hist(zscore,histbins);
hist(zscore,b);
hold on
%linepos = threshold*stdsumcounts+meansumcounts;
line([threshold threshold],ylim,'Color','r')
xlabel('Z-score Sum Spike Counts Across Channels')
ylabel('Num Bins')
subplot(3,1,2)
hist((maxcount-meansumcounts)/stdsumcounts,b)
hold on
%linepos = threshold*stdsumcounts+meansumcounts;
line([threshold threshold],ylim,'Color','r')
xlabel('Max Z-score Sum Counts per Trial')
ylabel('Num Bins')
subplot(3,1,3)
scatter(zscore(zscore<threshold),trialnums(zscore<threshold),'b');
hold on
scatter(zscore(zscore>=threshold),trialnums(zscore>=threshold),'r');
ylabel('Trial Number')
xlabel('Z-score Sum Spike Counts Across Channels')




trialstoremove = unique(trialnums(((sumcounts-meansumcounts)/stdsumcounts)>threshold));
bcidat(trialstoremove) = [];
end
