function [success, msgStr, fixWinOutput] = joystickMove(~, ~, positionHold, distanceTolerance, pixBoxLimit, cursorObjId, cursorR, cursorColor, immediateFailOnReturn)
% This function checks if the joystick has moved from a specified hold position. 

% success (1) if the joystick has moved from the hold position, noted by positionHold
% failure (0) if the joystick has not moved from the hold position

% cursor color
cursorColorDisp = [cursorColor(1) cursorColor(2) cursorColor(3)];

% compute the new cursor position and draw it on the display screen
[xVal, yVal, ~, ~] = sampleHallEffectJoystick();
xVal = xVal*pixBoxLimit;
yVal = yVal*pixBoxLimit;
cursorPos = [xVal, yVal]; 
cursorPosDisp = round(cursorPos); % round to prevent display computer from erroring
msgStr = sprintf('set %i oval 0 %i %i %i %i %i %i', [cursorObjId cursorPosDisp(1) cursorPosDisp(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);

% check how far cursor is from hold position
positionXHold = positionHold(1);
positionYHold = positionHold(2);
distanceFromHoldLoc = sqrt(sum(([xVal, yVal] - [positionXHold, positionYHold]).^2));
xyLargeEnough = checkWithinTolerance(distanceFromHoldLoc, 0, distanceTolerance, false);

if xyLargeEnough
    success = 1;
else
    success = 0;
end

% display on control computer 
cursorRadDisp = 5;
yellow  = [255 255 0];
winColors = yellow;

fixWinOutput = {[positionXHold cursorPosDisp(1)], [positionYHold cursorPosDisp(2)], [distanceTolerance cursorRadDisp], winColors};
