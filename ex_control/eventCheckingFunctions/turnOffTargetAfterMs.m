function [success, msgStr] = turnOffTargetAfterMs(loopStart,loopNow, timeBeforeTargOff, targetObjId)
% success if the cursor reaches the target
% failure if the cursor *doesn't* reach the target
global codes
persistent targetOff

targOff = false;
timeInLoopMs = 1000*(loopNow - loopStart);
if timeInLoopMs > timeBeforeTargOff
    targOff = true;
else
    targetOff = false;
end

% turn off target if it's too close
% NOTE AS OF RIGHT NOW THIS FUNCTION CAN'T TURN THE TARGET BACK ON!
if targOff
    targOffRadius = 0;
    targX = 0; % this isn't necessary per se, as the radius shuts off the target anyway
    targY = 0; % this isn't necessary per se, as the radius shuts off the target anyway
    randColor = [0 0 0];
    msgStr = sprintf('set %d oval 0 %i %i %i %i %i %i', targetObjId, [targX, targY, targOffRadius, randColor(1), randColor(2), randColor(3)]);
    if ~targetOff
        sendCode(codes.TARG_OFF)
        disp('TARGET OFF')
    end
    targetOff = true;
else
    msgStr = '';
end

% this is a nonblocking function for the loop--it's always successful!
success = 1;