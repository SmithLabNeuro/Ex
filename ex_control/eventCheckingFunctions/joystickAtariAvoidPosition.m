function success = joystickAtariAvoidPosition(~, ~, joystickButtonX, joystickButtonY, joystickDevInd)
% success if the joystick is not in the specified position
% failure -- none, only abort below
% abort if the joystick reaches specified position

joystickXY = Gamepad('GetAxis', joystickDevInd, 3:4);

% check for position
if isnan(joystickButtonX) && isnan(joystickButtonY)
    % allows us to ignore this joystick
    success = 1;
elseif all(joystickXY == [joystickButtonX, joystickButtonY])
    % joystick is in bad position
    success = -1;
else
    % joystick is not in bad position
    success = 1;
end

