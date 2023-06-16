function [success, msgStr, fixWinOutput] = bciShrinkAnnulus(loopStart, loopNow, bciSockets, e)
% success if the annulus size reduces to 0

global codes
persistent distToTargetState annulusRad holdAnnulusStart

if isempty(annulusRad)
    % At start of trial, no time bins are considered correct
    annulusRad = e.maxAnnulusRad;
    distToTargetState = 1;
end

if isempty(holdAnnulusStart)
    holdAnnulusStart = loopNow; % Set it to current time if empty
end

% for the control computer to track where target is
purple  = [255 0 255];
winColors = purple;

% annulus color (shown to subject on display computer)
col = e.annulusColor;
annulusColorDisp = [col(1) col(2) col(3)];

receivedMsg = '';
% this while loop ensures we only get the latest update from the BCI
% computer, as opposed to being stuck in previous ones, which might cause
% jumpiness if the control computer isn't fast enough
while matlabUDP2('check', bciSockets.sender)
    receivedMsg = matlabUDP2('receive', bciSockets.sender);
end
% Check if received annulus size from BCI computer 
if ~isempty(receivedMsg) && ~strcmp(receivedMsg, 'ack')
    try
        receivedMsgFromBci = typecast(uint8(receivedMsg), 'double')';
    catch err
        b = err;
        keyboard
    end
    % BCI computer sends us neural distance and annulus
    distToTargetState = receivedMsgFromBci(11);
    % Set annulusRad to received annulus radius
    annulusRad = receivedMsgFromBci(12);
end

% Need to send the BCI annulus radius at each time point
sendCode(codes.BCI_CURSOR_POS);
disp(annulusRad)
sendCode(annulusRad+5000); % ranges from 5000-6000
%sendCode(round(distToTargetState*1000));

success = 0;
% Check timing difference between holdAnnulusStart and loopNow
loopDiffMs = 1000*(loopNow - holdAnnulusStart); % Check that the annulus is maintained for given period of time
% Compare current annulus radius to a threshold to decide if success or fail: 
if (annulusRad < e.tolerance) && (loopDiffMs > e.msToMaintain)
    % Trial is only successful if past the threshold for a certain period
    % of time
    success = 1;
else
    % Reset timer if annulus is not smaller than tolerance
    if ~(annulusRad < e.tolerance)
        holdAnnulusStart =  loopNow;
    end
end

% send out the draw annulus commands whenever this function is called--is
% especially important for the fixation windows otherwise other functions
% will erase these windows when the cursor doesn't get updated
annX = e.centerX;
annY = e.centerY;
msgStr = sprintf('set %d annulus 0 %i %i %i %i %i %i %i',[e.objID, annX annY annulusRad e.thickness annulusColorDisp(1) annulusColorDisp(2) annulusColorDisp(3)]);

fixWinOutput = {[annX annX], [annY annY], [e.tolerance annulusRad], winColors};