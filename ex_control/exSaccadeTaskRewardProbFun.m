function [clicks,thisbehav] = exSaccadeTaskRewardProbFun(thisbehav,e)

% things to fix
% cannot do variable length deterministic intervals accurately



%% fields
if ~isfield(thisbehav,'isJackpot')
    thisbehav(1).isJackpot = [];
end

jackpots = [thisbehav.isJackpot];
%% define jackpot probabilities
if length(e.highJackpotProb)>1
    if length(jackpots)<=length(e.highJackpotProb)
        thishighjackpotprob = e.highJackpotProb(length(jackpots)+1);
    else
        thishighjackpotprob = e.highJackpotProb(end);
    end
else
    thishighjackpotprob = e.highJackpotProb;
end
 %       display(['Jackpot prob lengths: ',num2str(length(e.jackpotProb))])

if length(e.jackpotProb)>1
    if (length(jackpots)+1)<=length(e.jackpotProb)
%        display(['Jackpot Length: ',num2str(length(jackpots))])
        thisjackpotprob = e.jackpotProb(length(jackpots)+1);
    else
        thisjackpotprob = e.jackpotProb(end);
%        display(['Jackpot Length: ',num2str(length(jackpots))])
    end
else
    thisjackpotprob = e.jackpotProb;
end

%% define interval
if ~isfield(thisbehav,'jackpotInterval')
    thisbehav(1).jackpotInterval = [];
end

% the following allows for variable intervals
allintervals = [thisbehav.jackpotInterval];
if length(thisbehav)>1 && ~isempty(allintervals)
    thisbehav(end).jackpotInterval = allintervals(end);
else
    thisbehav(end).jackpotInterval = e.jackpotInterval;
end


%% decide if is jackpot
thisinterval = thisbehav(end).jackpotInterval;
%display(['jackpot: ',num2str(jackpots)])
% display(length(jackpots)>=thisinterval )
% display(thisinterval>1)
% if length(jackpots)>=thisinterval
% display(sum(jackpots((end-(thisinterval-1)+1):end))==0)
% end
%display(['This jackpot prob: ',num2str(thisjackpotprob)])
if length(jackpots)>=(thisinterval-1) && thisinterval>1 && mod((length(jackpots)+1),thisinterval)==0%sum(jackpots((end-(thisinterval-1)+1):end))==0   
    % skip jackpot randomly if previous jackpot trial was not a skip
    if isfield(e,'skipprob') && length(jackpots)>(2*thisinterval-1)&& jackpots(end-thisinterval+1)==1&&rand(1)<e.skipprob 
        thisbehav(end).isJackpot = 0;
    else
        thisbehav(end).isJackpot = rand(1)<=thisjackpotprob;
    end
elseif thisinterval==1
    thisrand = rand(1);
 %   display(['This Rand: ',num2str(thisrand)])
    thisbehav(end).isJackpot = thisrand <=thisjackpotprob;
 %   display(['isjackpot: ',num2str(thisbehav(end).isJackpot)])
else
    thisbehav(end).isJackpot = rand(1)<=thishighjackpotprob;
end

%% define upcoming interval
if thisbehav(end).isJackpot == 1
    thisbehav(end).jackpotInterval = e.jackpotInterval;
end



%% define reward amount
if length(e.jackpotClicks)>1
    if length(jackpots)<=length(e.jackpotClicks)
        thisjackpotclick = e.jackpotClicks(length(jackpots)+1);
    else
        thisjackpotclick = e.jackpotClicks(end);
    end
else
    thisjackpotclick = e.jackpotClicks;
end

if length(e.nonJackpotClicks)>1
    if length(jackpots)<=length(e.nonJackpotClicks)
        thisnonjackpotclick = e.nonJackpotClicks(length(jackpots)+1);
    else
        thisnonjackpotclick = e.nonJackpotClicks(end);
    end
else
    thisnonjackpotclick = e.nonJackpotClicks;
end




if thisbehav(end).isJackpot == 1
    clicks = thisjackpotclick;
else
    clicks = thisnonjackpotclick;
end
thisbehav(end).clicks = clicks;

end