function dat = removeBadTrials(dat)
baddat = zeros(length(dat),1);
for n = 1:length(dat)
    if length(dat(n).result)>1
        baddat(n) = 1;
        warning(['Removing dat from ' num2str(dat(n).time(1)),' to ',num2str(dat(n).time(2)),' because too many result codes'])
    end
end
dat = dat(baddat~=1);
end