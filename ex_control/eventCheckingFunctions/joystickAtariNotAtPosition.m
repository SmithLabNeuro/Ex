function [success] = joystickAtariNotAtPosition(~, ~, joystickButtonX, joystickButtonY, joystickDevInd)
% success if the joystick is not at specified position
% failure if the joystick gets into specified position
% abort -- none

joystickXY = Gamepad('GetAxis', joystickDevInd, 3:4);

% check for position
if isnan(joystickButtonX) && isnan(joystickButtonY)
    % allows us to ignore this joystick
    success = 1;
elseif all(joystickXY == [joystickButtonX, joystickButtonY])
    % joystick is at bad position
    success = 0;
else
    % joystick is not at bad position
    success = 1;
end

