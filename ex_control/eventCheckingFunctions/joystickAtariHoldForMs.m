function success = joystickAtariHoldForMs(loopStart, loopNow, msHold, joystickButtonX, joystickButtonY, joystickDevInd)
% success if the joystick *doesn't* change position for given time
% failure -- none, only abort below
% abort if the joystick *changes* position within given time

joystickXY = Gamepad('GetAxis', joystickDevInd, 3:4);

% check for a position change
if any(joystickXY - [joystickButtonX, joystickButtonY])
    success = -1;
else
    loopDiffMs = 1000*(loopNow-loopStart);
    if loopDiffMs > msHold
        success = 1;
    else
        success = 0;
    end
end

