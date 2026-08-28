function [success, msgStr, fixWinOutput] = joystickHold(loopStart, loopNow, positionXHold, positionYHold, distanceTolerance, pixBoxLimit, cursorObjectId, cursorR, cursorColorDisp, drawCursor, msHold)
% This function checks if the joystick is held in the correct position for the specified time. 

% success (1) if the joystick holds the correct position for the correct time, noted by positionXHold and positionYHold
% wait (0) if the joystick remains in position but the hold time has not yet passed
% failure (-1) if the joystick does not hold the correct position

success = 0;
[xVal, yVal, ~, ~] = sampleHallEffectJoystick();
xVal = xVal * pixBoxLimit;
yVal = yVal * pixBoxLimit;

% check that x, y, and z position are at desired hold positions
distanceFromHoldLoc = sqrt(sum(([xVal, yVal] - [positionXHold, positionYHold]).^2));
xySmallEnough = checkWithinTolerance(distanceFromHoldLoc, 0, distanceTolerance, true);

% if the hold time has passed, success no matter what (so this function
% doesn't care what you do past the hold time)
loopDiffMs = 1000*(loopNow-loopStart);
if loopDiffMs > msHold
    success = 1;
else
    if ~xySmallEnough
        success = -1; % failure if the joystick does not hold the correct position
    else
        success =  0; % wait if the joystick remains in position but the hold time has not yet passed
    end
end

% now we want to draw on both the control computer and the display computer
yellow  = [255 255 0];
winColors = yellow;
cursorPos = [xVal, yVal]; 

% this communicates with showex to update the cursor position on the display computer
cursorPosDisp = round(cursorPos); % round to prevent display computer from erroring
if drawCursor
    % draw the cursor    
    msgStr = sprintf('set %d oval 0 %i %i %i %i %i %i', [cursorObjectId cursorPosDisp(1) cursorPosDisp(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);
else
    msgStr = '';
end

% this draws on the control computer to show where the cursor is relative to the hold position
numWindows = 2;
maxSizeInfoVals = 2;
sizeInfo = nan(maxSizeInfoVals, numWindows);
sizeInfo(1:length(distanceTolerance),1) = distanceTolerance;
sizeInfo(1:length(cursorR),2) = cursorR;
fixWinOutput = {[positionXHold cursorPosDisp(1)], [positionYHold cursorPosDisp(2)], sizeInfo,winColors};

end