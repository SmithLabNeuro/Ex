function [fixSuccess, freqSuccess, freqvals,timevec ]= waitForMSmatlabUDP(waitTime,fixX,fixY,targX,targY,r,targetfreq, freqtolerance,freqperflag,freqperthresh,freqdurflag,freqdurfun,cursorvisibleflag,varargin)
% function success = waitForMS(waitTime,fixX,fixY,r)
% 
% ex trial helper function: waits for t ms, checking to ensure that the eye
% remains within the fixation window.  If time expires, fixSuccess = 1, 
% but if fixation is broken first, fixSuccess returns 0
%
% waitTime: time to maintain fixation (in ms)
% fixX, fixY: in pixels, the offset of the fixation from (0,0)
% r: in pixels, the radius of the fixation window
%
% 2015/08/14 by Adam Snyder and Matt Smith. Now allows user to pass a
% 'recenterFlag' such that the fixX and fixY will be ignored and instead
% the current eye position is used.
%
% 2016/5/24 Ryan Williamson added functionality to check on success based
% on digital codes

global params; 
global codes;
global sockets;
global bciCursorTraj; 
timevec = [];

checkbciflag = 0;
if ~exist('freqperflag')
    freqperflag = 0;
elseif freqperflag == 1
    if ~exist('freqperthresh')
        error('Frequency percentage threshold undefined when flag is on')
    end
end
 % check inputs
    winColors = [255 255 0];
    bciloopinterval = 0.04996;
    cursoroffsettime = 0.02;
    recenterFlag = false;
    
    if nargin > 4        
        vx = 1;
        while vx <= numel(varargin),
            switch class(varargin{vx})
                case 'char'
                    recenterFlag = varargin{vx+1};      
                    vx = vx+2;
                otherwise
                    if ~isempty(varargin{vx}),
                        winColors = varargin{vx};
                    end;
                    vx = vx+1;
            end;
        end;   
    end;
    
    if recenterFlag,
        d=samp;
        eyePos = projectCalibration(d(end,:)); %changed to new project calibration function (supports polynomial regressors) -ACS 29Oct2013
        fixX = eyePos(1);
        fixY = eyePos(2); % removed fixX/Y - MAS Sept 2019
    end;
    
    if nargin >= 4
        if params.bciCursorEnabled
            bciColors = [0 255 0];
            drawFixationWindows([fixX targX targetfreq(1)],[fixY targY targetfreq(2)],[r r freqtolerance],[winColors;winColors;bciColors]);
        else
            drawFixationWindows([fixX],[fixY],[r]',[winColors]);
        end
    elseif nargin ~=1
        error('waitForMS can have exactly 1, 4 or 5 input arguments');
    end
    
    fixSuccess = 1;
    
%     while matlabUDP2('check',sockets(2))
%          matlabUDP2('receive',sockets(2));
%     end
%     
    if cursorvisibleflag == 1
        cursorsize = 50;
       % msg('set 4 oval 0 %i %i %i %i %i %i',[0 0 cursorsize 0 255 0]);
    %    msg('set 5 oval 0 %i %i %i %i %i %i',[0 0 cursorsize 0 255 0]);
       
    %    msg('obj_switch 4 -5')
        msgAndWait(['diode 4']);   
        msgAndWait(['obj_on 4']); 
        %msg('obj_on 4')
        objon = 4;
        sendCode(codes.STIM1_ON);
    end
    thisStart = tic;
    
    if nargin >= 4 
      
         freqSuccess = 0;
         successclockflag = 0;
         freqrewardflag = 0;
         freqvals = [0 0];
         
        bciupdatetime = GetSecs;
        cursorupdatetime = bciupdatetime + cursoroffsettime;
        bciind =0;
        while (toc(thisStart)*1000) <= waitTime
            loopTop = GetSecs;
            d=samp;
           
            eyePos = projectCalibration(d(end,:)); %changed to new project calibration function (supports polynomial regressors) -ACS 29Oct2013

 %           eyePos = eyePos - [fixX fixY];

            relPos = bsxfun(@minus,eyePos(:),[fixX;fixY]); %position relative to each window        
            switch size(r,1)
                case 1, %circular window        
                    inWin = sum(relPos.^2,1)<r.^2; 
                case 2, %rectangular window
                    inWin = all(abs(relPos)<abs(r),1);
                otherwise
                    error('EX:waitForMS:badRadius','Radius must have exactly 1 or 2 rows');
            end;
            % retrieve digital events from Ripple and check freq
           freqevents = [];
            if GetSecs - bciupdatetime >= bciloopinterval
                checkbciflag = 1;
                bciupdatetime = GetSecs;
            end
%             if bciind == 5
%                 if matlabUDP2('check',sockets(2))
%                     freqvals = [freqvals; str2double(strsplit(matlabUDP2('receive',sockets(2))))];
%                 end
%                 bciind = 0;
%             else
%                 bciind = bciind + 1;
%             end
%             
%             msgAndWait('dset 4 oval 0 %i %i %i %i %i %i',[freqvals(end,1) freqvals(end,2) cursorsize 0 255 0]);
%             if bciind == 0
%                 timevec = [timevec; (toc(thisStart)*1000)];
%             end
%             sendCode(codes.STIM1_ON)
%             
            if checkbciflag == 1
                if matlabUDP2('check',sockets(2))
                    %bciupdatetime = GetSecs;
                    while matlabUDP2('check',sockets(2))
                        freqevents = [freqevents; str2double(strsplit(matlabUDP2('receive',sockets(2))))];
                        %timevec = [timevec; (toc(thisStart)*1000)];
                        freqvals = [freqvals; freqevents];
                        sendCode(codes.SOUND_CHANGE);
                    end 
                    checkbciflag = 0;
                end
            end
            if GetSecs - cursorupdatetime >= bciloopinterval
                
                if ~isempty(freqvals)
                 cursorupdatetime = bciupdatetime + cursoroffsettime;   
                    % send info to display computer
                    
                    bciCursorTraj = freqvals;
                    
                    if cursorvisibleflag == 1
                        sendCode(codes.STIM1_ON)
                        msgAndWait('dset 4 oval 0 %i %i %i %i %i %i',[freqvals(end,1) freqvals(end,2) cursorsize 0 255 0]);
                        sendCode(codes.STIM2_ON)
                        %                     if objon == 4
                            
%                             msg('set 5 oval 0 %i %i %i %i %i %i',[freqvals(end,1) freqvals(end,2) cursorsize 0 255 0]);
%                             msg('obj_switch -4 5')
%                             sendCode(codes.STIM2_ON)
% %                             msg('obj_on 5')
% %                             msg('obj_off 4')
%                             objon = 5;
%                         else
%                             msg('set 4 oval 0 %i %i %i %i %i %i',[freqvals(end,1) freqvals(end,2) cursorsize 0 255 0]);
%                             msg('obj_switch 4 -5')
% %                             msg('obj_on 4')
% %                             msg('obj_off 5')
%                             objon = 4;
%                             sendCode(codes.STIM1_ON)
%                         end
                    end
                    timevec = [timevec; (toc(thisStart)*1000)];
                end
                
            end
%             if ~isempty(freqevents)
%                     % save freqs in behavioral data
%                     
%                     freqvals = [freqvals; freqevents];
%                     
%                     
%                     bciCursorTraj = freqvals;
                  
                    
                    % check if freq is within threshold of target. For now this
                    % does not require the target to be held 
                   
%                     if ~isempty(find(tolerancefunction(freqevents,targetfreq,freqtolerance),1))%~isempty(find(abs(freqevents-targetfreq)<freqtolerance,1))
%                         if successclockflag == 0
%                             successclock =  1;
%                             successclockflag = 1;
%                         elseif successclockflag==1 && successclock < targettime
%                             successclock = successclock + 1;
%                         end
%                         if successclockflag && successclock >= targettime
%                             freqSuccess = 1;
%                             successclockflag = 0;
%                             if(juicedur>0)
%                                 giveJuice(1,juicedur,juicedur);
%                                 sendCode(codes.REWARD);
%                             end
%                         end
%                     else
%                         successclockflag = 0;
%                     end
               
%             end
            
            if freqdurflag == 1              
                if ~isempty(freqvals) && freqdurfun(freqvals,targetfreq)==1
                    fixsuccess = 1;
                    freqsuccess = 1;
                    %break
                end
            end
            
            if keyboardEvents()||~inWin
                fixSuccess = 0;
                break;
            end
            if (GetSecs-loopTop)>params.waitForTolerance, warning('waitFor:tooSlow','waitForMS exceeded latency tolerance - %s',datestr(now)); end; %warn tolerance exceeded -acs22dec2012
            
        end
    else %don't worry about fixation window - this is essentially just a pause (can be broken with a key press)
        while (toc(thisStart)*1000) <= waitTime
            loopTop = GetSecs;            
            if keyboardEvents()
                fixSuccess = 0;
                break;
            end
            if (GetSecs-loopTop)>params.waitForTolerance, warning('waitFor:tooSlow','waitForMS exceeded latency tolerance - %s',datestr(now)); end; %warn tolerance exceeded -acs22dec2012
        end %waitTime
    end
    
    if nargin > 2
        drawFixationWindows()
    end   
    if freqperflag == 1
        if (sum(freqvals>targetfreq)/length(freqvals))>freqperthresh
            freqSuccess = 1;
        end
    end
end
