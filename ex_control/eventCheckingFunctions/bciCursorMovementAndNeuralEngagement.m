function [success, msgStr, fixWinOutput, extraFuncOutput] = bciCursorMovementAndNeuralEngagement(~, loopNow, bciSockets, binSizeMs, initCursorPos, initVelocity, initNeuralEngagement, pixBoxLimit, targX,targY,targRadius, cursorObjId, cursorSizeInfo, cursorColor, targWinCursRad, timestepsInTargetForCorrect)
% success if the cursor reaches the target
% failure if the cursor *doesn't* reach the target

global codes
persistent cursorPos velocityCurr loopTimeLast neCurr timestepsCorrect 

center = [0; 0];
vecCenterToTarget = [targX-center(1); targY - center(2)];
normVecCenterToTarget = vecCenterToTarget/norm(vecCenterToTarget);
normVecOrth = null(normVecCenterToTarget');

% centerOutAssistVel = normVecCenterToTarget*centerToOutVelForAssist;
if isempty(cursorPos)
    cursorPos = initCursorPos;
    velocityCurr = initVelocity;  
    neCurr = initNeuralEngagement;
    loopTimeLast = loopNow;
    timestepsCorrect = 0;
end

% for the control computer to track where target is
purple  = [255 0 255];
winColors = purple;
binChange = false;

% cursor color
cursorColorDisp = [cursorColor(1) cursorColor(2) cursorColor(3)];

% cursor size
cursorRectHalfWidth = cursorSizeInfo(1);
cursorRectHalfHeight = cursorSizeInfo(2);

% these are the 0 values around which we'll determine cursor position (for
% saving to the NEV file and also the MAT file)
posX0 = 10000;
posY0 = 10000;
posNE0 = 10000;
% posShiftForCode = [posX0 posY0 posEng];

% grab the loop time delta for the move
msPerS = 1000;
binSizeS = binSizeMs/msPerS;
loopTimeDiff = min(loopNow - loopTimeLast, binSizeS);
% disp(loopTimeDiff)

receivedMsg = '';
% this while loop ensures we only get the latest update from the BCI
% computer, as opposed to being stuck in previous ones, which might cause
% jumpiness if the control computer isn't fast enough
ind = 0;
while matlabUDP2('check', bciSockets.sender)
    receivedMsg = matlabUDP2('receive', bciSockets.sender);
    ind = ind+1;
end
if ind>1
disp(ind)
disp(ind)
disp(ind)
disp(ind)
disp(ind)
end
if ~isempty(receivedMsg) && ~strcmp(receivedMsg, 'ack')
    % compute the new cursor position (and don't let it get out of the bounding
    % box)
    try
        velocityAndNeuralEngagementNew = typecast(uint8(receivedMsg), 'double')';
        velocityFromBci = velocityAndNeuralEngagementNew(1:2);
        neCurr = velocityAndNeuralEngagementNew(3);
        binChange = true;
        %             disp(neCurr)
    catch err
        b = err;
        keyboard
    end
    % constrain the velocity based on automonkey/passive assist
    centerToTargetVel = velocityFromBci' * normVecCenterToTarget;
    orthVel = velocityFromBci' * normVecOrth;
    velocityCurr = centerToTargetVel*normVecCenterToTarget + orthVel*normVecOrth;
    
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

% the below should rest at 45 degrees (when neNew=0), and one standard
% deviation above and below (since neNew is z-scored) will span -54 to 54
% degrees; basically, it can give good feedback for +/- 1 z-scored NE
% change--the stretchFactor will allow us to increase what z-score span
% moves between those angles; 
% as written, +/- 2 z-score to span most of 0-90 degrees (if this wasn't there, it
% would go between ~10-80 degrees, where 0 and 90 get reached at +/-
% infinity
stretchFactor = 1;
cursorAngleDeg = 1.2*atand(neCurr/stretchFactor)/2;
% making sure the cursor angle remains between -45 and 45 degrees
cursorAngleDeg = min(cursorAngleDeg, 45);
cursorAngleDeg = max(cursorAngleDeg, -45);
            
% send the constrained velocity back to the BCI computer; no heads up
% messages or receive accept messages so that the loop can be tight--but
% the BCI computer better be ready to receive what we're sending!
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

posShiftX = posX0 + cursorPos(1);
posShiftY = posY0 + cursorPos(2);
neShift = posNE0 + sign(neCurr)*min(abs(neCurr)*100, 500);
sendCode(codes.BCI_CURSOR_POS);
sendCode(posShiftX);
sendCode(posShiftY);
sendCode(neShift);

% compute how close the cursor is to the target
relPos = (cursorPos-[targX; targY]);
switch size(targRadius,1)
    case 1 %circular window
        distToTarget = sqrt(sum(relPos.^2));
        targetHit = distToTarget < targWinCursRad;
    case 2 %rectangular window
        targetHit = all(abs(relPos)<abs(targRadius+targWinCursRad),1);
    case 3 % rotated rectangle
%         % the point here is to ask whether the center of mass of the target
%         % is within the rotated target rectangle--easiest way I could think
%         % of was to unrotate the target rectangle (and the relevant
%         % relative position), and do the normal rectangle check
%         targUnrotAngle = -targRadius(3);
%         unrotMatTarg = [cosd(targUnrotAngle), -sind(targUnrotAngle);
%                   sind(targUnrotAngle), cosd(targUnrotAngle)];
%                
%         relCursPtsTargRef = unrotMatTarg * relPos;
%         success = all(abs(relCursPtsTargRef)<abs(targRadius(1:2)+targWinCursRad),1);
        % the point here is for the cursor to be completely swallowed by
        % the target in order for him to get reward.
        targUnrotAngle = -targRadius(3);
        unrotMatTarg = [cosd(targUnrotAngle), -sind(targUnrotAngle);
                  sind(targUnrotAngle), cosd(targUnrotAngle)];
        rotMatCursor = [cosd(cursorAngleDeg), -sind(cursorAngleDeg);
                  sind(cursorAngleDeg), cosd(cursorAngleDeg)];
               
        relCursPtsTargRef = unrotMatTarg * relPos;
        rectanglePtsRelTarg = relCursPtsTargRef + ...
            unrotMatTarg*rotMatCursor*[-cursorRectHalfWidth  cursorRectHalfHeight;
              cursorRectHalfWidth  cursorRectHalfHeight;
              cursorRectHalfWidth -cursorRectHalfHeight;
             -cursorRectHalfWidth -cursorRectHalfHeight]'; % transpose, don't forget
        minX = min(rectanglePtsRelTarg(1, :));
        maxX = max(rectanglePtsRelTarg(1, :));
        minY = min(rectanglePtsRelTarg(2, :));
        maxY = max(rectanglePtsRelTarg(2, :));
        
        targetHit = minX>-targRadius(1) && maxX < targRadius(1) &&...
            minY>-targRadius(2) && maxY<targRadius(2);

    otherwise
        error('EX:waitForFixation:badRadius','Radius must have exactly 1 or 2 rows');
end

if targetHit && binChange
    % add an in-target timepoint when he's in the target
    disp('change bin!!')
    disp('change bin!!')
    timestepsCorrect = timestepsCorrect+1;
else
    % reset timesteps to zero if he ever leaves target
    timestepsCorrect = 0;
end

success = targetHit && (timestepsCorrect >= timestepsInTargetForCorrect);
loopTimeLast = loopNow;

% send out the draw cursor commands whenever this function is called--is
% especially important for the fixation windows otherwise other functions
% will erase these windows when the cursor doesn't get updated
cursorPosDisp = round(cursorPos); % round to prevent display computer from erroring
cursorAngleDisp = round(cursorAngleDeg);
msgStr = sprintf('set %i rotpoly 0 %i %i %i %i %i %i %i %i', [cursorObjId, cursorPosDisp(1) cursorPosDisp(2) cursorRectHalfWidth cursorRectHalfHeight cursorAngleDisp cursorColorDisp(1) cursorColorDisp(2) cursorColorDisp(3)]);
% msgStr
fixWinOutput = {[targX cursorPosDisp(1)], [targY cursorPosDisp(2)], [targWinCursRad cursorRectHalfWidth], winColors};
extraFuncOutput.cursorPosFinal = cursorPos;
extraFuncOutput.velocityFinal = velocityCurr;
extraFuncOutput.cursorAngleFinal = cursorAngleDeg;



