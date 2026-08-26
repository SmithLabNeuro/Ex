function [success, msgStr, fixWinOutput] = joystickMoveToTarget(~,~, targX, targY, cursorObjectId, cursorR, cursorColor, targWinRad, pixBoxLimit)
% This function allows cursor movement to a target location. 
% It checks if the cursor is within a specified radius of the target location.
% It also updates cursor position on the display computer and provides information for the control computer to track the cursor and target locations.

% success (1) if the cursor is within the target radius, noted by targX and targY
% failure (0) if the cursor is not within the target radius

% for the control computer to track where target is
yellow  = [255 255 0];
winColors = yellow;

% cursor color
cursorColorDisp = [cursorColor(1) cursorColor(2) cursorColor(3)];

% compute the new cursor position (and don't let it get out of the bounding box)
[xVal, yVal, ~, ~] = sampleHallEffectJoystick();
cursorPos = [xVal*pixBoxLimit, yVal*pixBoxLimit]; 
signCursor = sign(cursorPos);
cursorPos = signCursor.*min(pixBoxLimit, abs(cursorPos));

% compute how close the cursor is to the target
relPos = ([targX targY] - cursorPos);
distToTarget = sqrt(sum(relPos.^2));

success = distToTarget < targWinRad; % 1 or 0

% draw the cursor
cursorPosDisp = round(cursorPos); % round to prevent display computer from erroring
msgStr = sprintf('set %d oval 0 %i %i %i %i %i %i', [cursorObjectId cursorPosDisp(1) cursorPosDisp(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);

fixWinOutput = {[targX cursorPosDisp(1)], [targY cursorPosDisp(2)], [targWinCursRad cursorR], winColors};
