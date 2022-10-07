function rewardamount = slowdriftRewardFun(driftsmooth,largepupilrewardflag)

numdevs = sign(driftsmooth)*floor(abs(driftsmooth));
if largepupilrewardflag == 1
    numdevs = abs(max(0,numdevs))+1;   
else
    numdevs = abs(min(0,numdevs))+1;
end

points = 1:numdevs;
rewarddist = [];
for n = 1:length(points)
    rewarddist = [rewarddist points(n)*ones(1,numdevs-n+1)];
end

rewardamount = rewarddist(randi(length(rewarddist)));
end