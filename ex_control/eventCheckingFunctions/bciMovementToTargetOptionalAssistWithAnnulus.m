function [success, msgStr, fixWinOutput, extraFuncOutput] = bciMovementToTargetOptionalAssistWithAnnulus(~, loopNow, bciSockets, binSizeMs, initCursorPos, initVelocity, pixBoxLimit, targX,targY,targRadius, cursorObjId, cursorR, cursorColor, targWinCursRad, centerToTargetBciScale, orthogonalBciScale, centerToOutVelForAssist, e)
% success if the cursor reaches the target
% failure if the cursor *doesn't* reach the target

global codes
persistent cursorPos velocityCurr loopTimeLast numBin

maxAnnulusSize = e.maxAnnulusRad;
annulusColorDisp = e.annulusColor;
cursorOpacity=e.cursorOpacity;
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
    numBin = 0;
    cursorPos = initCursorPos;
    velocityCurr = initVelocity;
    if all(velocityCurr == 0)
        if any(abs(centerOutAssistVel) > 0)
            velocityCurr = centerOutAssistVel;
        else
            % this forces the BCI computer to start at zero
            % velocity--useful for starting after a freeze period
            try
                uint8Msg = typecast(velocityCurr, 'uint8');
            catch err
                a = err;
                keyboard
            end
        end
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


receivedMsg = '';
% this while loop ensures we only get the latest update from the BCI
% computer, as opposed to being stuck in previous ones, which might cause
% jumpiness if the control computer isn't fast enough
ind2 = 0;
while matlabUDP2('check', bciSockets.sender)
    receivedMsg = matlabUDP2('receive', bciSockets.sender);
    ind2 = ind2+1;
end
if ind2>1
end
if ~isempty(receivedMsg) && ~strcmp(receivedMsg, 'ack')
    % compute the new cursor position (and don't let it get out of the bounding
    % box)
    try
        velocityFromBci = typecast(uint8(receivedMsg), 'double')';
        numBin = numBin + 1;
    catch err
        b = err;
        keyboard
    end
    % constrain the velocity based on automonkey/passive assist
    centerToTargetVel = centerToTargetBciScale * velocityFromBci' * normVecCenterToTarget;
    orthVel = orthogonalBciScale * velocityFromBci' * normVecOrth;
    velocityCurr = centerToTargetVel*normVecCenterToTarget + orthVel*normVecOrth + centerOutAssistVel;
    
else
    velocityFromBci = velocityCurr;
end


posPixelChange = velocityCurr * loopTimeDiff;
cursorPosNew = cursorPos + posPixelChange;
signCursor = sign(cursorPosNew);
cursorPosLimited = signCursor.*min(pixBoxLimit, abs(cursorPosNew));
posPixelChangeLimit = cursorPosLimited - cursorPos;
if ~all(abs(posPixelChangeLimit - posPixelChange)<1e-10)
    velocityCurr = posPixelChangeLimit/loopTimeDiff;
    cursorPos = cursorPosLimited;
else
    cursorPos = cursorPosNew;
end

% send the constrained velocity (from both passive assist and border
% limitations back to the BCI computer; no heads up messages or receive
% accept messages so that the loop can be tight--but the BCI computer
% better be ready to receive what we're sending!
if any(~ismembertol(velocityFromBci, velocityCurr))
    try
        uint8Msg = typecast(velocityCurr, 'uint8');
    catch err
        a = err;
        keyboard
    end
    sendMessageWaitAck(bciSockets, uint8Msg');
end



% I think this posShift, a double, gets cast to an int in
% unixSendByte, which sendCode calls... keep in mind in case weird
% things appear, but it's been working until now...
posShift = posShiftForCode + cursorPos;
posShiftX = posShift(1);
posShiftY = posShift(2);
sendCode(codes.BCI_CURSOR_POS);
sendCode(posShiftX);
sendCode(posShiftY);

 % compute how close the cursor is to the target
 relPos = ([targX; targY] - cursorPos);
 
 maxTargDistance = sqrt(sum([targX; targY].^2));
 switch size(targRadius,1)
     case 1 %circular window
         distToTarget = sqrt(sum(relPos.^2));
         success = distToTarget < targWinCursRad;
     case 2 %rectangular window
         success = all(abs(relPos)<abs(targRadius),1);
     case 3 % rotatedcursorOpacity=100;
 rectangle
         success = all(abs(relPos)<abs(targRadius),1);
     otherwise
         error('EX:waitForFixation:badRadius','Radius must have exactly 1 or 2 rows');
 end


loopTimeLast = loopNow;

% send out the draw cursor commands whenever this function is called--is
% especially important for the fixation windows otherwise other functions
% will erase these windows when the cursor doesn't get updated
cursorPosDisp = round(cursorPos); % round to prevent display computer from erroring


msgStrOne = sprintf('set %i oval 0 %i %i %i %i %i %i %i', [cursorObjId, cursorPosDisp(1) cursorPosDisp(2) cursorR cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3) cursorOpacity]);

fixWinOutput = {[targX cursorPosDisp(1)], [targY cursorPosDisp(2)], [targWinCursRad cursorR], winColors};

extraFuncOutput.cursorPosFinal = cursorPos;
extraFuncOutput.velocityFinal = velocityCurr;

% Annulus size should be the same as the cursor distance

disp([distToTarget, maxTargDistance, maxAnnulusSize])
distRatio = (distToTarget - targWinCursRad)/(maxAnnulusSize - e.minAnnulusRad);
if distRatio <= 0 
    annulusRad = e.minAnnulusRad;
elseif distRatio >= 1
    annulusRad = maxAnnulusSize;
else
    % Annulus needs to be an integer value or else showex will error
    annulusRad = round(distRatio*maxAnnulusSize);    
end

% send out the draw annulus commands whenever this function is called--is
% especially important for the fixation windows otherwise other functions
% will erase these windows when the cursor doesn't get updated
annX = e.centerX;
annY = e.centerY;
msgStrTwo = sprintf('set %d annulus 0 %i %i %i %i %i %i %i',[e.objID, annX annY annulusRad e.thickness annulusColorDisp(1) annulusColorDisp(2) annulusColorDisp(3)]);
msgStr = [msgStrOne, ' ', msgStrTwo];
%disp(msgStr)