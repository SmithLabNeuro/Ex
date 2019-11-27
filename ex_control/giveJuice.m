function giveJuice(juiceX,juiceInterval,juiceTTLDuration,randomizeReward)
% function giveJuice(juiceX,juiceInterval,juiceTTLDuration)
%
% provides a reward.  uses juiceX and juiceTTLDuration to determine the
% value of reward. If no arguments are provided for either one, uses the
% global params values.
%
% juiceX - number of juice rewards
% juiceInterval - time in ms between rewards
% juiceTTLDuration - if on computer control, ms length of juice reward (set
%     to 1 ms if not on computer control)
%
% Modified 19Feb2013 by Adam Snyder to enable randomization of parameters
% (and send a structure with the parameters to the data machine and the
% allCodes file).
%
% Modified 4Dec2015 by Matt Smith to make it flexible to run unixSendPulse
% (on linux) or winReward (on a PC) to deliver the juice
%
% Modified 11Aug2017 by Matt Smith to remove windows capability

global params;

if params.rewarding
    %Adding hard-coded limits on the parameters here, since we can now draw the
    %parameters from a distribution with potentially infinite support.
    limitsJuiceTTLDuration  = [1 500];  %in ms
    limitsJuiceX            = [1 10];   %number of clicks
    limitsJuiceInterval     = [1 500]; %in ms
    %Setting defaults:
    if (nargin < 4), randomizeReward = false; end; %the randomizeReward
    if (nargin < 3), juiceTTLDuration = params.juiceTTLDuration; end;
    if (nargin < 2), juiceInterval = params.juiceInterval; end;
    if (nargin < 1), juiceX = params.juiceX; end;
    if randomizeReward||isfield(params,'randomizeReward'),
        if randomizeReward||params.randomizeReward,
            if isfield(params,'rewardRandomizationParameter'), %valid values are 'juiceTTLDuration', 'juiceX', or 'juiceInterval'
                rewardRandomizationParameter = params.rewardRandomizationParameter;
            else
                rewardRandomizationParameter = 'juiceTTLDuration'; %default
            end;
            if isfield(params,'rewardDistribution'),
                rewardDistribution = params.rewardDistribution;
            else
                rewardDistribution = 'Normal'; %default
            end;
            if isfield(params,'rewardDistributionParams') 
                rewardDistributionParams = num2cell(params.rewardDistributionParams);
            else
                val = eval(rewardRandomizationParameter);
                rewardDistributionParams = {val, 2*sqrt(val)};
            end;
            randVal = round(random(rewardDistribution,rewardDistributionParams{:})); %#ok<NASGU> %restrict to integers
            eval(sprintf('%s = randVal',rewardRandomizationParameter));
            checkLimits;
            sendStruct(struct('juiceX',juiceX,'juiceTTLDuration',juiceTTLDuration,'juiceInterval',juiceInterval));
        end;
    else
        checkLimits;
    end;    
    
    % new comedi-based Linux digital output function
    unixSendPulse(params.juiceChan,juiceTTLDuration);
    
    if isfield(params,'rewardSound')&&params.rewardSound,
        if params.soundEnabled
            playTone(500,100);
        end
        % old code from windows sound commands
        %         t = 0:0.0001:0.05;      % 10 seconds @ 1kHz sample rate
        %         fo = 200; f1 = 400;   % Start at 10Hz, go up to 400Hz
        %         y = chirp(t,fo,max(t),f1,'logarithmic');
        %         sound(repmat(0.5*y,1,3),10000);
    end;
    
    for i = 2:juiceX
        pause(juiceInterval/1000); % wait for the next reward
        % new comedi-based Linux digital output function
        unixSendPulse(params.juiceChan,juiceTTLDuration);
    end
end


    function checkLimits
        juiceTTLDuration(juiceTTLDuration<min(limitsJuiceTTLDuration))=min(limitsJuiceTTLDuration); juiceTTLDuration(juiceTTLDuration>max(limitsJuiceTTLDuration))=max(limitsJuiceTTLDuration);
        juiceX(juiceX<min(limitsJuiceX))=min(limitsJuiceX); juiceX(juiceX>max(limitsJuiceX))=max(limitsJuiceX);
        juiceInterval(juiceInterval<min(limitsJuiceInterval))=min(limitsJuiceInterval); juiceInterval(juiceInterval>max(limitsJuiceInterval))=max(limitsJuiceInterval);
    end
        
end
