function [success, msgStr, fixWinOutput] = joystickHallEffectGrabCheck(~, ~, positionXHold, positionYHold, angleZHold, distanceTolerance, angleTolerance, pixelDistForMaxJoystickPos)
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
% angLargeEnough = zValAng > (angleZHold-angleTolerance);
angLargeEnough = abs(zValAng) > (angleZHold-angleTolerance);%checkWithinTolerance(zValAng, angleZHold, -angleTolerance, false);

% if button
%     disp('button held! button held! button held!')
%     rewardTimeProportionHold = 0.1; % reward for just holding...
%     juiceX = 1; % number of rewards
%     juiceInterval = 1; % time between rewards (ignored if juiceX = 1)
%     juiceTTLDuration = round(rewardTimeProportionHold*params.juiceTTLDuration); % duration of reward
%     giveJuice(juiceX,juiceInterval,juiceTTLDuration);
%     
% end

prStr = sprintf('angle %6.2f\n', zValAng);
prStrB = prStr(1:end);
prStrB(:) = sprintf('\b');
fprintf(prStrB)
fprintf(prStr)

if (xySmallEnough && angLargeEnough)%buttonPress || 
%     button = true;
    success = 1;
end

cursorR = 5;
cursorIndicatorAngle = 20; % degree wedge
cursAngDispWedgeVals = [cursorR, cursorIndicatorAngle/2, zValAng]';

yellow  = [255 255 0];
winColors = yellow;
msgStr = '';
numWindows = 3;
maxSizeInfoVals = 3;
sizeInfo = nan(maxSizeInfoVals, numWindows);
sizeInfo(1:length(distanceTolerance),1) = distanceTolerance;
sizeInfo(1:length(cursorR),2) = cursorR;
sizeInfo(1:length(cursAngDispWedgeVals),3) = cursAngDispWedgeVals;

% drawFixationWindows([positionXHold cursorPosDisp(1)], [positionYHold cursorPosDisp(2)], [distanceTolerance cursorR],winColors);
% fixWinOutput = {[positionXHold cursorPosDisp(1)], [positionYHold cursorPosDisp(2)], [distanceTolerance cursorR],winColors};
fixWinOutput = {[positionXHold cursorPosDisp(1) cursorPosDisp(1)], [positionYHold cursorPosDisp(2) cursorPosDisp(2)], sizeInfo,winColors};

