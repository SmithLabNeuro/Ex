function success = joystickAtariInteract(~, ~, joystickNoInteractX, joystickNoInteractY, joystickDevInd, checkAll)
% success if the *joystick* (not the cursor) reaches correct position
% failure if the joystick does not reach correct position
if nargin < 6
    % assumes checkAll is not provided and sets default value
    checkAll = true;
end

successCheck = zeros(length(joystickDevInd), 1);
for i=1:length(joystickDevInd)
    currDevInd = joystickDevInd(i);
    joystickXY = Gamepad('GetAxis', currDevInd, 3:4);
    % check for position being correct
    if sum(abs(joystickXY - [joystickNoInteractX, joystickNoInteractY]))
        successCheck(i) = 1;
    else
        successCheck(i) = 0;
    end
end
if checkAll
    success = all(successCheck);
else
    success = any(successCheck);
end
