function success = joystickAtariOrHallEffectInteract(loopStart, loopTop, atariInputs, heInputs)
% this function checks for either one or the other joystick being clicked

% this global is to keep a history of the used joystick
global joystickUsedHistory

successAtari = joystickAtariInteract(loopStart, loopTop, atariInputs{:});
successHE = joystickHallEffectInteract(loopStart, loopTop, heInputs{:});

success = successAtari | successHE;
if success
    joystickUsedHistory = [joystickUsedHistory; [successAtari, successHE]];
end

