function success = joystickTemplate(loopStart, loopNow, cursX, cursY, cursR, joystickDevInd, cursorColor, cursorR, joystickPxPerSec)

persistent posCurr

if isempty(posCurr)
    posCurr = [cursX, cursY];
end


loopDiff = loopStart - loopNow;
joystickMaxValue = 32768;
adjFactorForVel = joystickPxPerSec*loopDiff/joystickMaxValue;

posMv = Gamepad('GetAxis', joystickDevInd, 3:4);
posCurr = posCurr + posMv*adjFactorForVel;

if any(posMv ~= 0)
    success = -1;
end

