function [success, msgStr, fixWinOutput, extraOut] = joystickHallEffectUntwistCheck(~, ~, initAngle, positionXHold, positionYHold, angleZHold, distanceTolerance, angleTolerance, pixelDistForMaxJoystickPos)
% success if the *joystick* (not the cursor) reaches correct position
% failure if the joystick does not reach correct position

success = 0;
[xVal, yVal, zValAng, buttonPress] = sampleHallEffectJoystick();
% if ~buttonPress
%     disp('huzzah?')
% end
pixBoxLimit = pixelDistForMaxJoystickPos;
xVal = xVal * pixBoxLimit;
yVal = yVal * pixBoxLimit;

cursorPos = [xVal, yVal]; 
cursorPosDisp = round(cursorPos); % round to prevent display computer from erroring

% global button

% check that x, y, and z position are at desired hold positions
distanceFromHoldLoc = sqrt(sum((cursorPos - [positionXHold, positionYHold]).^2));
xySmallEnough = checkWithinTolerance(distanceFromHoldLoc, 0, distanceTolerance, true);
% the check for untwisting is that it has passed some required untwist
% stage *or* gone past it--the first half of the logical checks whether
% it's within the untwist stage, but maybe not past it while the second
% half also is correct if it's past it (i.e. went from a positive to a
% negative twist or vice versa)
angSmallEnough = (abs(zValAng) < (angleZHold+angleTolerance)) || (sign(zValAng) ~= sign(initAngle));%checkWithinTolerance(zValAng, angleZHold, -angleTolerance, false);


prStr = sprintf('angle %6.2f\n', zValAng);
prStrB = prStr(1:end);
prStrB(:) = sprintf('\b');
fprintf(prStrB)
fprintf(prStr)

if (xySmallEnough && angSmallEnough)
    success = 1;
end

cursorR = 5;
cursorIndicatorAngle = 20; % degree wedge
yellow  = [255 255 0];
winColors = yellow;
msgStr = '';
cursAngDispWedgeVals = [cursorR, cursorIndicatorAngle/2, zValAng]';

numWindows = 3;
maxSizeInfoVals = 3;
sizeInfo = nan(maxSizeInfoVals, numWindows);
sizeInfo(1:length(distanceTolerance),1) = distanceTolerance;
sizeInfo(1:length(cursorR),2) = cursorR;
sizeInfo(1:length(cursAngDispWedgeVals),3) = cursAngDispWedgeVals;

fixWinOutput = {[positionXHold cursorPosDisp(1) cursorPosDisp(1)], [positionYHold cursorPosDisp(2) cursorPosDisp(2)], sizeInfo, winColors};
extraOut = {zValAng}; % so that we can do something with how much he twists away
