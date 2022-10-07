function [success, msgStr] = turnOffTargetWhenCursorClose(~,~, targX,targY,targRadius, targOffWinCursRad, pixelDistForMaxJoystickPos, targetObjId)
% success if the cursor reaches the target
% failure if the cursor *doesn't* reach the target



pixBoxLimit = pixelDistForMaxJoystickPos;


% compute the new cursor position (and don't let it get out of the bounding
% box)
[xVal, yVal, ~, ~] = sampleHallEffectJoystick();
cursorPos = [xVal*pixBoxLimit, yVal*pixBoxLimit]; 

% compute how close the cursor is to the target
relPos = ([targX targY] - cursorPos);
distToTarget = sqrt(sum(relPos.^2));

switch size(targRadius,1)
    case 1 %circular window
        turnOff = distToTarget < targOffWinCursRad;
        
%         if targY==0 && targX > 0 && distToTarget < targWinCursRad
%             disp(distToTarget)
%         end
    case 2 %rectangular window
        turnOff = all(abs(distToTarget)<abs(targRadius),1);
    otherwise
        error('EX:waitForFixation:badRadius','Radius must have exactly 1 or 2 rows');
end


% turn off target if it's too close
% NOTE AS OF RIGHT NOW THIS FUNCTION CAN'T TURN THE TARGET BACK ON!
if turnOff
    targOffRadius = 0;
    randColor = [0 0 0];
    msgStr = sprintf('set %d oval 0 %i %i %i %i %i %i', targetObjId, [targX, targY, targOffRadius, randColor(1), randColor(2), randColor(3)]);
else
    msgStr = '';
end

% this is a nonblocking function for the loop--it's always successful!
success = 1;

