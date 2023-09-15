function [success, msgStr, fixWinOutput] = autoShrinkAnnulus(loopStart, loopNow, e)
% success if the annulus size reduces to 0

global codes
persistent annulusRad holdAnnulusStart 
if isempty(annulusRad)
    annulusRad = e.maxAnnulusRad;
end

if isempty(holdAnnulusStart)
    % Set both to current time if empty
    holdAnnulusStart = loopStart; 
end

% for the control computer to track where target is
purple  = [255 0 255];
winColors = purple;

% annulus color (shown to subject on display computer)
col = e.annulusColor;
annulusColorDisp = [col(1) col(2) col(3)];

success = 0;
% Check timing difference between holdAnnulusStart and loopNow
loopDiffMs = 1000*(loopNow - holdAnnulusStart); % Check that the annulus is maintained for given period of time
% Interpolating the annulus values so that it gradually shrinks at a
% constant rate for entire delay period.
numSteps=round(e.delayMs/e.autoMonkeyAnnulusBinSize);
interpAnnulusRad=linspace(e.maxAnnulusRad,e.tolerance,numSteps);
interpAnnulusIdx=round(loopDiffMs/e.autoMonkeyAnnulusBinSize);
if(interpAnnulusIdx<1)
    interpAnnulusIdx=1;
end
if(interpAnnulusIdx>=numSteps)
    interpAnnulusIdx=numSteps;
end
annulusRadSend= round(interpAnnulusRad(interpAnnulusIdx));
if(interpAnnulusIdx>=numSteps)
    success=1
end
% Send the Annulus Cursor position
sendCode(codes.BCI_CURSOR_POS);
sendCode(annulusRad+5000); % ranges from 5000-6000
% send out the draw annulus commands whenever this function is called--is
% especially important for the fixation windows otherwise other functions
% will erase these windows when the cursor doesn't get updated
annX = e.centerX;
annY = e.centerY;
msgStr = sprintf('set %d annulus 0 %i %i %i %i %i %i %i',[e.objID, annX annY annulusRadSend e.thickness annulusColorDisp(1) annulusColorDisp(2) annulusColorDisp(3)]);
fixWinOutput = {[annX annX], [annY annY], [e.tolerance annulusRad], winColors};
