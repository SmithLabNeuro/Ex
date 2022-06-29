function rewardamount = slowdriftRewardFun2(driftsmooth,largepupilrewardflag,rewardpercentage)
if ~exist('rewardpercentage','var')
    rewardpercentage = 0.2;
end


numdevs = sign(driftsmooth)*floor(abs(driftsmooth));
if largepupilrewardflag == 1
    numdevs = abs(max(0,numdevs))+1;   
else
    numdevs = abs(min(0,numdevs))+1;
end

randnum = rand(1)<rewardpercentage;
if randnum
    rewardamount = numdevs;
else
    rewardamount = 1;
end
end