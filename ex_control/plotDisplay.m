function plotDisplay(obj,~)
%function plotDisplay
%
% called by a timer initiated by runex. plots the current eye position on
% the screen.  
%
% Modified:
%
% 2017/08/15 by Adam Snyder: changed from using Screen('CopyWindow') to
% Screen('DrawTexture').
%
% 2015/12/22 by Matt Smith: Now uses samp for eye/mouse position, and
% buffers the last position in "lastSamp" to allow it to plot the eye
% position line. Also removed the "mid" calculation and replaced with a
% transpose. Not sure why that was there.
%
% 2013/09/12 by Adam Snyder: removed references to histogram window.

global wins bciCursorTraj params;
persistent lastSamp;

if ~exist('lastSamp','var')
    lastSamp = samp;
end

% buffer the last eye position so we can draw a line
newSamp = samp;
eyeVolts = [lastSamp;newSamp]; 
lastSamp = newSamp;

% convert voltages to pixel space (relative to the whole screen)
eyePix = projectCalibration(eyeVolts);

% invert the y-value for display purposes
eyePix(:,2) = -eyePix(:,2);

% convert to the right scale for the pixel display and the voltage display
dispPix = bsxfun(@plus,eyePix * diag(wins.pixelsPerPixel),wins.midE);
dispVolts = bsxfun(@plus,eyeVolts * diag(wins.pixelsPerVolt),wins.midV);

% left eye position figure - voltage space
Screen('DrawLines',wins.voltage,dispVolts',1,wins.eyeTraceColor);
Screen('DrawTexture',wins.w,wins.voltage,[0 0 wins.voltageSize(3:4)-wins.voltageSize(1:2)],wins.voltageSize);
Screen('FillOval',wins.w,wins.eyePosColor,[dispVolts(end,:) - wins.eyeDotRad, dispVolts(end,:) + wins.eyeDotRad]+[wins.voltageSize(1:2) wins.voltageSize(1:2)]);




% right eye position figure - pixel space
Screen('DrawLines',wins.eye,dispPix',1,wins.eyeTraceColor);
Screen('DrawTexture',wins.w,wins.eye,[0 0 wins.eyeSize(3:4)-wins.eyeSize(1:2)],wins.eyeSize);
if params.bciCursorEnabled && params.bciEnabled && ~isempty(bciCursorTraj)
    % Ryan stuff here
    tempbciCursorTraj = bciCursorTraj;
    tempbciCursorTraj(:,2) = -1*tempbciCursorTraj(:,2);
    scaledbciCursorTraj =bsxfun(@plus,tempbciCursorTraj * diag(wins.pixelsPerPixel),wins.midE);
    Screen('DrawLines',wins.eye,scaledbciCursorTraj',1,[0 255 0]);
    %Screen('DrawTexture',wins.w,wins.eye,[0 0 wins.eyeSize(3:4)-wins.eyeSize(1:2)],wins.eyeSize);
    Screen('FillOval',wins.w,[0 255 0],[scaledbciCursorTraj(end,:) - 2*wins.eyeDotRad, scaledbciCursorTraj(end,:) + 2*wins.eyeDotRad]+[wins.eyeSize(1:2) wins.eyeSize(1:2)]);
end
Screen('FillOval',wins.w,wins.eyePosColor,[dispPix(end,:) - wins.eyeDotRad, dispPix(end,:) + wins.eyeDotRad]+[wins.eyeSize(1:2) wins.eyeSize(1:2)]);


    
Screen('DrawTexture',wins.w,wins.info,[0 0 wins.infoSize(3:4)-wins.infoSize(1:2)],wins.infoSize);

% IMPORTANT: need to use the '2' below in the Screen('Flip') call to make sure 
% plotDisplay doesn't hold up other runex activities.
%
% Modified by MAS in December 2015 to use the '2' for the "dontsync" argument
% and this solved reaction time quantization. This was part of the porting
% of "Ex" to linux. Also need to make sure the updating of this plotter
% isn't too fast. 30 ms should be good.
Screen('Flip',wins.w,0,0,2); 

end
