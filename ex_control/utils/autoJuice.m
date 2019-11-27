function autoJuice(interval,nRewards)
% function autoJuice(interval,nRewards)
%  interval: in seconds
%  nRewards: total number of rewards (all single juice)

for i = 1:nRewards
    if ispc
        % old Windows-compatible code
        winReward;
    else
        % new comedi-based Linux digital output function
        unixSendPulse(params.juiceChan);
    end
    pause(interval); % wait 150 ms before the 2nd reward
end
