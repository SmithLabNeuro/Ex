function success = cursorHold(~, ~, joystickDevInd)
% success if the cursor *doesn't* change position
% failure -- none, only abort below
% abort if the cursor *changes* position (really, this is failure if the
% joystick is anywhere but in the center)

posMv = Gamepad('GetAxis', joystickDevInd, 3:4);

if any(posMv ~= 0)
    success = -1;
else
    success = 1;
end

