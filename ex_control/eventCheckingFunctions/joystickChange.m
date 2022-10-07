function success = joystickChange(~, ~, joystickButtonX, joystickButtonY, joystickDevInd)
% success if the joystick *changes* position
% failure if the joystick *doesn't* change position

% variable so that once success if achieved, it stays achieved throughout
% the waitForEvent
persistent successAchieved

if ~isempty(successAchieved) && successAchieved
    % once successful always successful
    success = 1;
else
    % just check if the joystick has changed position at all
    joystickXY = Gamepad('GetAxis', joystickDevInd, 3:4);
    
    if ~any(joystickXY - [joystickButtonX, joystickButtonY])
        success = 0;
    else
        success = 1;
        % ensure that success remains even if the joystick goes back to the
        % original position
        successAchieved = 1;
    end
end