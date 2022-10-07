function success = joystickHallEffectHoldCheck(~, ~, positionXHold, positionYHold, angleZHold, distanceTolerance, angleTolerance, pixelDistForMaxJoystickPos)
% success if the *joystick* (not the cursor) stays at correct position/twist
% failure if the joystick does not stay correct position/twist

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

if ~(xySmallEnough && angLargeEnough)
    success = -1;
else
    success = 1;
end

end

function passTolCheck = checkWithinTolerance(actual, expected, tolerance, checkInIfTrue)
    passTolCheck = false;
    % flag to see if the one must stay within the tolerance or be outside
    % the tolerance
    checkOutIfTrue = ~checkInIfTrue;
    if checkInIfTrue && (abs(actual - expected) < tolerance)
        passTolCheck = true;
    elseif checkOutIfTrue &&  (abs(actual-expected) > tolerance)
        passTolCheck = true;
    end
end