function success = joystickAtariHold(~, ~, joystickButtonX, joystickButtonY, joystickDevInd)
% success if the joystick *doesn't* change position
% failure -- none, only abort below
% abort if the joystick *changes* position

joystickXY = Gamepad('GetAxis', joystickDevInd, 3:4);

% check for a position change
if any(joystickXY - [joystickButtonX, joystickButtonY])
    success = -1;
else
    success = 1;
end

