function rewardamount = slowdriftRewardFun2(driftsmooth,largepupilrewardflag,rewardpercentage)
if ~exist('rewardpercentage','var')
    rewardpercentage = 0.2;
end


numdevs = sign(driftsmooth)*floor(abs(driftsmooth));
if largepupilrewardflag == 1
    numdevs = abs(max(0,numdevs))+2;   %the + constant sets the minimum value of jackpot trials
else
    numdevs = abs(min(0,numdevs))+2;
end

randnum = rand(1)<rewardpercentage;
if randnum
    rewardamount = numdevs;
else
    rewardamount = 1; % this constant sets the value of non-jackpot trials
end
end