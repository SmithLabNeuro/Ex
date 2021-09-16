function success = joystickAtariAcquirePositionOrHallEffectGrabCheck(loopStart, loopTop, atariInputs, heInputs)
% this function checks for either one or the other joystick being clicked

% this global is to keep a history of the used joystick
global joystickUsedHistory

successAtari = joystickAtariAcquirePosition(loopStart, loopTop, atariInputs{:});
successHE = joystickHallEffectGrabCheck(loopStart, loopTop, heInputs{:});

success = successAtari | successHE;
if success
    joystickUsedHistory = [joystickUsedHistory; [successAtari, successHE]];
end

