function success = joystickGetClick(~, ~, analogChJoystickButton)

joystickButtonState = unixGetAnalogInput(analogChJoystickButton);

success = joystickButtonState;