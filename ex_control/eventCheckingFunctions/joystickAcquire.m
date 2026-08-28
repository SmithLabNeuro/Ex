function [success, msgStr, fixWinOutput] = joystickAcquire(~, ~, positionXHold, positionYHold, distanceTolerance, pixBoxLimit)
% This function checks if the joystick is in the correct position to acquire a target. 
% It does not draw anything on the display computer, and it does not update any object locations. 
% It simply checks if the joystick is within a certain distance of a specified hold position.

% success (1) if the joystick reaches correct position, noted by positionXHold and positionYHold
% failure (0) if the joystick does not reach correct position

success = 0;
[xVal, yVal, ~, ~] = sampleHallEffectJoystick(); % returns in volts

% convert to pixels
xVal = xVal * pixBoxLimit;
yVal = yVal * pixBoxLimit;

% check that x, y position are at desired hold positions
distanceFromHoldLoc = sqrt(sum(([xVal, yVal] - [positionXHold, positionYHold]).^2));
xySmallEnough = checkWithinTolerance(distanceFromHoldLoc, 0, distanceTolerance, true);

if xySmallEnough
    success = 1;
end

% unused outputs, not drawing or updating object locations
msgStr = '';
fixWinOutput = {};