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

global wins bciCursorTraj params typingNotes notes sqlDb;
persistent lastSamp lastClick updateNotes;

if isempty(updateNotes)
    updateNotes = true;
end

if ~exist('lastSamp','var')
    lastSamp = samp;
end
if isempty(lastClick)
    lastClick = tic;
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


% info screen
Screen('DrawTexture',wins.w,wins.info,[0 0 wins.infoSize(3:4)-wins.infoSize(1:2)],wins.infoSize);

% % lab notebook screen
if ~isempty(sqlDb)
    if isfield(params, 'keyboardName')
        clickTextExpansion = 1.5;
        keyboardId = GetKeyboardIndices([params.keyboardName]);
        [pos(1), pos(2), mouseClick] = GetMouse();
        if mouseClick(1)
            timeSinceLastClick = toc(lastClick);
            if timeSinceLastClick > 0.3 % at least half a second between clicks...
                lastClick = tic;
                pos = (pos-(wins.labNotesSize(1:2) + wins.labNotesDim(3:4)./2)) .*[1 1];
                posInWindowBottomSixth = pos(2) > wins.labNotesDim(4)/3 && pos(2) < wins.labNotesDim(4)/2 && pos(1) > -wins.labNotesDim(3)/2 && pos(1) < wins.labNotesDim(3)/2;% should be bottom sixth of screen
                
                if posInWindowBottomSixth
                    if ~typingNotes
                        disp('start')
                        for ind = 1:length(keyboardId)
                            KbQueueCreate(keyboardId(ind)); % creates a queue for the keyboard
                            KbQueueStart(keyboardId(ind)); % starts the queue for the keyboard
                            disp(ind)
                        end
                        typingNotes = true;
                    else
                        disp('end')
                        typingNotes = false;
                        updateNotes = true;
                        % update the session with any notes that have been written
                        sqlDb.exec(sprintf('UPDATE experiment_session SET notes = "%s" WHERE session_number = %d AND animal = "%s"', notes, params.sessionNumber, params.SubjectID));
                        % NOTE: you might think the below is a good idea, but
                        % it is not. It kills the ability to do any inputs to
                        % the trial. Kthx.
                        % for ind = 1:length(keyboardId)
                        %     KbQueueStop(keyboardId(ind));
                        % end
                    end
                end
            end
        end
        if ~typingNotes && updateNotes
            updateNotes = false;
            % start typing button
            typeNotesButtonText = 'Click here to type notes';
            x=0;y=wins.labNotesDim(4)-clickTextExpansion*wins.textSize;
            gray=(WhiteIndex(0)+BlackIndex(0))/2;
            Screen(wins.labNotes,'FillRect',gray);
            Screen('TextSize',wins.labNotes,clickTextExpansion*wins.textSize);
            Screen('DrawText',wins.labNotes,typeNotesButtonText,x,y, [], gray);
            
            % fill in what notes are there...
            x=0;y=wins.labNotesDim(4)-(1+clickTextExpansion)*wins.textSize;
            Screen('TextSize',wins.labNotes,wins.textSize);
            wrappedString = WrapString(char(notes), 35);
            drawStr = split(wrappedString, newline);
            for lns = length(drawStr):-1:1
                Screen('DrawText',wins.labNotes,drawStr{lns},x,y, [], gray);
                y = y - wins.textSize;
            end
        elseif typingNotes
            % prep the background
            white=WhiteIndex(0);
            Screen(wins.labNotes,'FillRect',white);
            
            % stop typing button
            stopTypeNotesButtonText = 'Click here to stop typing';
            x=0;y=wins.labNotesDim(4)-clickTextExpansion*wins.textSize;
            gray=(WhiteIndex(0)+BlackIndex(0))/2;
            Screen('TextSize',wins.labNotes,clickTextExpansion*wins.textSize);
            Screen('DrawText',wins.labNotes,stopTypeNotesButtonText,x,y, [], gray);
            
            notesInd = length(notes)+1;
            for ind = 1:length(keyboardId)
                [event,nremaining] = KbEventGet(keyboardId(ind));
                if ~isempty(event)
                    if event.CookedKey==sprintf('\b')
                        notesInd = notesInd-1;
                        if notesInd==0
                            notesInd = 1;
                        else
                            notes(notesInd) = [];
                        end
                    else
                        notes(notesInd) = event.CookedKey;
                        notesInd = notesInd+1;
                    end
                end
                while nremaining
                    [event,nremaining] = KbEventGet(keyboardId(ind));
                    if event.CookedKey==sprintf('\b')
                        notesInd = notesInd-1;
                        if notesInd==0
                            notesInd = 1;
                        else
                            disp('huh')
                            notes(notesInd) = [];
                        end
                    else
                        notes(notesInd) = event.CookedKey;
                        notesInd = notesInd + 1;
                    end
                end
                notes = notes(notes~=0);
                notes(notes==13) = newline;
                %             disp(notes)
                
                %             if ~isempty(notes)
                x=0;y=wins.labNotesDim(4)-(1+clickTextExpansion)*wins.textSize;
                Screen('TextSize',wins.labNotes,wins.textSize);
                wrappedString = WrapString(char(notes), 35);
                drawStr = split(wrappedString, newline);
                for lns = length(drawStr):-1:1
                    Screen('DrawText',wins.labNotes,drawStr{lns},x,y, [], white);
                    y = y - wins.textSize;
                end
                %             end
            end
        end
    end
    Screen('DrawTexture',wins.w,wins.labNotes,[0 0 wins.labNotesSize(3:4)-wins.labNotesSize(1:2)],wins.labNotesSize);
end

% IMPORTANT: need to use the '2' below in the Screen('Flip') call to make sure 
% plotDisplay doesn't hold up other runex activities.
%
% Modified by MAS in December 2015 to use the '2' for the "dontsync" argument
% and this solved reaction time quantization. This was part of the porting
% of "Ex" to linux. Also need to make sure the updating of this plotter
% isn't too fast. 30 ms should be good.
Screen('Flip',wins.w,0,0,2); 

end
