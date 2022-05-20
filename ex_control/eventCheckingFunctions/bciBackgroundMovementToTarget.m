function [success, msgStr, fixWinOutput, extraFuncOutput] = bciBackgroundMovementToTarget(~, loopNow, bciSockets, binSizeMs, initCursorPos, initVelocity, pixBoxLimit, targX,targY, cursorR, cursorColor, targWinCursRad, centerToTargetBciScale, orthogonalBciScale, centerToOutVelForAssist)
% success if the cursor reaches the target
% failure if the cursor *doesn't* reach the target

global codes
persistent cursorPos velocityCurr loopTimeLast


if nargin < 10
    centerToTargetBciScale = 1;
    orthogonalBciScale = 1;
    centerToOutVelForAssist = 0; % this is only used if the two scales are 0
end

center = [0; 0];
vecCenterToTarget = [targX-center(1); targY - center(2)];
normVecCenterToTarget = vecCenterToTarget/norm(vecCenterToTarget);
normVecOrth = null(normVecCenterToTarget');

centerOutAssistVel = normVecCenterToTarget*centerToOutVelForAssist;

if isempty(cursorPos)
    cursorPos = initCursorPos;
    velocityCurr = initVelocity;
    if all(velocityCurr == 0)
        velocityCurr = centerOutAssistVel;
    end
    loopTimeLast = loopNow;
end

% for the control computer to track where target is
purple  = [255 0 255];
winColors = purple;

% cursor color
cursorColorDisp = [cursorColor(1) cursorColor(2) cursorColor(3)];

% these are the 0 values around which we'll determine cursor position (for
% saving to the NEV file and also the MAT file)
posX0 = 10000;
posY0 = 10000;
posShiftForCode = [posX0 posY0];

% grab the loop time delta for the move
msPerS = 1000;
binSizeS = binSizeMs/msPerS;
loopTimeDiff = min(loopNow - loopTimeLast, binSizeS);
% disp(loopTimeDiff)


success = 1; % this function just plots in the background, so no fails!
if matlabUDP2('check', bciSockets.sender)
    if isequal(velocityCurr, [0;0])
        disp(velocityCurr)
    end
    receivedMsg = matlabUDP2('receive', bciSockets.sender);
    if ~isempty(receivedMsg) && ~strcmp(receivedMsg, 'ack')
        % compute the new cursor position (and don't let it get out of the bounding
        % box)
        try
            velocityNew = typecast(uint8(receivedMsg), 'double')';
        catch err
            b = err;
            keyboard
        end
        % constrain the velocity based on automonkey/passive assist
        centerToTargetVel = centerToTargetBciScale * velocityNew' * normVecCenterToTarget;
        orthVel = orthogonalBciScale * velocityNew' * normVecOrth;
        velocityCurr = centerToTargetVel*normVecCenterToTarget + orthVel*normVecOrth + centerOutAssistVel;
        
        % send the constrained velocity back to the BCI computer; no heads
        % up messages or receive accept messages so that the loop can be
        % tight--but the BCI computer better be ready to receive what we're
        % sending!
        try
            uint8Msg = typecast(velocityCurr, 'uint8');
        catch err
            a = err;
            keyboard
        end
        sendMessageWaitAck(bciSockets, uint8Msg');
        
    else
        disp('message waiting but not received?')
    end
end


posPixelChange = velocityCurr * loopTimeDiff;
cursorPos = cursorPos + posPixelChange;
cursorPos = min(pixBoxLimit, cursorPos);


% I think this posShift, a double, gets cast to an int in
% unixSendByte, which sendCode calls... keep in mind in case weird
% things appear, but it's been working until now...
posShift = posShiftForCode + cursorPos;
posShiftX = posShift(1);
posShiftY = posShift(2);
sendCode(codes.BCI_CURSOR_POS);
sendCode(posShiftX);
sendCode(posShiftY);

loopTimeLast = loopNow;
        
% send out the draw cursor commands whenever this function is called--is
% especially important for the fixation windows otherwise other functions
% will erase these windows when the cursor doesn't get updated
cursorPosDisp = round(cursorPos); % round to prevent display computer from erroring
%         disp(size(cursorPos))
%         disp(cursorPosDisp(1))
%         disp(cursorPosDisp(2))
msgStr = '';

extraFuncOutput.cursorPosFinal = cursorPos;
extraFuncOutput.velocityFinal = velocityCurr;

fixWinOutput = {[targX cursorPosDisp(1)], [targY cursorPosDisp(2)], [targWinCursRad cursorR], winColors};
