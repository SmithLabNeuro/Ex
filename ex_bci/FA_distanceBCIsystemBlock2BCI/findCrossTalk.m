function [removechans,finalcorrmat] = findCrossTalk(filename)
tic
%maxlag = 0;
nev=readNEV(filename);
lasttime = max(nev(:,3));
bins = 0:1/10000:(1/10000+lasttime);

spikenev = nev(nev(:,1)~=0,:);
channels = unique(spikenev(:,1));
spikecell= sparse(length(bins)-1,length(channels));
meanvec = zeros(length(channels),1);
varvec = zeros(length(channels),1);
for n = 1:length(channels)
    n
    inds = spikenev(:,1)==channels(n);
    spikecell(:,n) = (histcounts(spikenev(inds,3),bins))';
    meanvec(n) = sum(inds)/length(inds);
    %varvec(n) = (sum((spikecell(:,n)-meanvec(n)).^2)/length(spikecell(:,n)))^0.5;
end
toc
tic
sparsecell = sparse(spikecell);
meanmat = zeros(length(channels),length(channels));
for m = 1:length(meanvec)
    display(m)
    meanmat(m,:) = sum(meanvec(m)*sparsecell,1);
end
%c=repmat((meanvec)',size(sparsecell,1),1);
%meanmat = sparsecell'*c;%sparsecell'*repmat(meanvec',size(sparsecell,1),1)
meanmat2 = 2*meanmat;
meansqmat = meanvec*meanvec'*size(spikecell,1);
xcorrmat = sparsecell'*sparsecell;
corrmat = full(xcorrmat)-meanmat2+meansqmat;

stdvec = diag(corrmat).^0.5;
stdmat = stdvec*stdvec';
finalcorrmat = corrmat./stdmat;


[a,b]=ind2sub(size(corrmat),find(tril(finalcorrmat,-1)>0.1));
c=[a(a~=b),b(a~=b)];
removechans = [];
while ~isempty(c)
numperchan = zeros(length(channels),1);
for n = 1:length(channels)
    numperchan(n) = sum(sum(c==n));
end
[~,maxind] = max(numperchan);
thischan = channels(maxind);
removechans = [removechans thischan];
c = c(sum(c == maxind,2)==0,:);
end
removechans = sort(removechans);
end 