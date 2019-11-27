function pos=samp(arg)
%function pos=samp(arg)
%
% Samples Eye Movement Data, returning current X/Y position.
% Currently grabs from PCIM-DAS1602 card if params.eyeTrackerAnalog is
% true, and if it's false then it should try to use the EyeLink Toolbox
% (not currently implemented, but rough code is there). In the analog case,
% uses the Linux comedi-based Mex files (unixGetEyes). 
%
% Note: As of November 2015, the histogram functions of samp are no
% longer supported. Samp is only for eye movements now.
%
% Modified 11Aug2017 by Matt Smith - windows features removed, linux only
%
% -------------------------------------------------------------------------
% OLD WINDOWS COMMENTS HERE:
%
% recompile with this command:
% mex winSamp.c cbw32.lib
%
% if background sampling is not running, samp always simply starts the timer and returns nothing.
%
% samp() returns the current eye position (x,y)
% samp(n) where n is positive returns the last n eye positions in a [n x 2] matrix
%
% samp(-3) marks the start of a histogram run and returns the index in the buffer (index for debugging only)
% samp(-2) marks the histogram align point and returns the index in the buffer (index for debugging only)
% samp(-1) marks the histogram end point and returns the index in the buffer (index for debugging only)
% samp(0) returns [h a], where h is the histogram from start to end
%
% samp(-4) stops the background sampling and frees the memory.
%
% NOTE: the buffer is 10 seconds long, so if the total time recording a hist
% exceeds this there will be unexpected results.  Make the buffer longer if you
% anticipate >10s trials.  This is in samp.c, under the variable "samples".  It
% MUST be a multiple of 3.  Right now it is 30000, which is 10000 samples for 3
% channels.  So for 30 seconds you'd make it 90000.  I don't know how this affects
% performance.  It may not.

global params eyeHistory eyeHistoryCurrentPos wins;
persistent eyelinkStatus eyelinkEyeUsed;

defaultPosVal = [-32768 -32768]; pos = defaultPosVal;

% if you're in mouse mode, samp should return the mouse position
if params.getEyes == 0
    [pos(1), pos(2)]=GetMouse;
%     sca;
    % converts top-left origin mouse coordinates to voltage-equivalent 
    % (and flip the y-coordinate so negative is down)
    pos = (pos-(wins.voltageDim(3:4)./2)) .*wins.voltScale./(wins.voltageDim(3:4)/2) .*[1 -1];
%     pos = (pos-(wins.voltageDim(3:4)./2)) .* (wins.voltScale/500)./(wins.voltageDim(3:4)/2) .*[1 -1];
else
    if params.eyeTrackerAnalog
        if nargout == 1
            pos=unixGetEyes();
        else
            return;
        end
    else
        % insert new eyelink toolbox code here
        % NOTE: Could consider coordinating calibration routine with Eyelink's
        % calibration. Would need separate code for that
        
        % start eye tracking if it's not started, return current position if it is
        error('Ethernet eye tracking via Eyelink Toolbox not supported yet');
%         
%         if nargin == 0
%             if (isempty(eyelinkStatus) || eyelinkStatus ~= 0)
%                 eyelinkStatus = Eyelink('Initialize');
%                 if (eyelinkStatus ~= 0)
%                     error('Can not initialize eyelink');
%                 end
%                 el=EyelinkInitDefaults;
%                 eyelinkEyeUsed = Eyelink('EyeAvailable');
%                 if eyelinkEyeUsed == el.BINOCULAR
%                     error('Binocular tracking not supported yet.');
%                 end
%                 eyelinkEyeUsed = eyelinkEyeUsed + 1; % indexing starts at zero
%                 
%                 % make sure that we get gaze data from the Eyelink
%                 % is this needed????
%                 %            Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');
%                 Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA,GAZERES,HREF,PUPIL,STATUS,INPUT,HMARKER');
%                 
%                 % open file to record data to
%                 edfFile='extemp.edf';
%                 Eyelink('Openfile', edfFile);
%                 
%                 Eyelink('StartRecording');
%             end
%             if nargout == 1
%                 if Eyelink('EyeAvailable') ~= -1
%                     evt=Eyelink('NewestFloatSample');
%                     pos = [evt.gx(eyelinkEyeUsed) evt.gy(eyelinkEyeUsed)];
%                 else % value to indicate eyes aren't being tracked
%                     pos = defaultPosVal;
%                 end
%             end
%         else
%             if arg==-4 % stop sampling and free everything
%                 if eyelinkStatus==0
%                     try
%                         Eyelink('StopRecording');
%                         Eyelink('CloseFile');
%                         %Eyelink('SetOfflineMode');
%                     catch
%                         disp('Problem trying to stop recording');
%                     end
%                     % set it to empty so it will force reinitialization
%                     eyelinkStatus = [];
%                 end
%             else
%                 error(['Eyelink error: arg ',num2str(arg),' not supported']);
%             end
%             if nargout == 1
%                 pos = defaultPosVal;
%             end
%         end
        
        %     if nargout == 1,
        %         pos = 5.*pos./3276.8; %rescale from int to volts -acs02dec2015
        %     end
    end
    
end

% store eye position in a global history buffer
eyeHistory(eyeHistoryCurrentPos,:) = [pos(1) pos(2) GetSecs];
eyeHistoryCurrentPos = mod(eyeHistoryCurrentPos,size(eyeHistory,1)) + 1;
if ~any(isnan(eyeHistory(:)))&&range(eyeHistory(:,3)*1000)<=params.eyeSmoothing,
    warning('samp:fullBuffer','The eye history buffer filled faster than the requested smoothing duration - increase buffer or shorten smoothing');
end;

% if smoothing is on, grab the position samples < X ms ago and average them
if params.eyeSmoothing > 1
    pos = nanmean(eyeHistory(((GetSecs-eyeHistory(:,3)) * 1000) < params.eyeSmoothing,1:2),1);
end

%if params.getEyes,
%    pos = pos .* wins.pixelsPerMV + wins.midV;
%    pos(2) = wins.voltageDim(4) - pos(2); % flip Y coordinate
%end;

%disp(nanmean(eyeHistory(:,1:2)))
%disp(max(eyeHistory(:,3))-min(eyeHistory(:,3)));

end
