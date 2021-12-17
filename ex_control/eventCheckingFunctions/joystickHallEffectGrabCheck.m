function success = joystickHallEffectGrabCheck(~, ~, positionXHold, positionYHold, angleZHold, distanceTolerance, angleTolerance, pixelDistForMaxJoystickPos)
% success if the *joystick* (not the cursor) reaches correct position
% failure if the joystick does not reach correct position

success = 0;
[xVal, yVal, zValAng, ~] = sampleHallEffectJoystick();
pixBoxLimit = pixelDistForMaxJoystickPos;
xVal = xVal * pixBoxLimit;
yVal = yVal * pixBoxLimit;

% check that x, y, and z position are at desired hold positions
distanceFromHoldLoc = sqrt(sum(([xVal, yVal] - [positionXHold, positionYHold]).^2));
xySmallEnough = checkWithinTolerance(distanceFromHoldLoc, 0, distanceTolerance, true);
% angLargeEnough = zValAng > (angleZHold-angleTolerance);
angLargeEnough = abs(zValAng) > (angleZHold-angleTolerance);%checkWithinTolerance(zValAng, angleZHold, -angleTolerance, false);

if xySmallEnough && angLargeEnough
    prStr = sprintf('angle %3.1f\n', zValAng);
    prStrB = prStr(1:end);
    prStrB(:) = sprintf('\b');
    fprintf(prStrB)
    fprintf(prStr)
    success = 1;
end
