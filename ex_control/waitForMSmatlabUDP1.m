function [fixSuccess, bciSuccess, sentcursors,timevec ]= waitForMSmatlabUDP1(waitTime,fixX,fixY,targX,targY,r,bcitarget, bcitolerance,bcifun,updatesintargetwin,cursorvisibleflag,cursorcolor,cursorsize,bcirewardflag,varargin)
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

 % check inputs
    winColors = [255 255 0];
    bciloopintervalpre = 0.04996;
    cursoroffsettime = 0.020;
    bciloopinterval = bciloopintervalpre - cursoroffsettime;
    recenterFlag = false;
    sentcursors = [];
    
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
            drawFixationWindows([fixX targX bcitarget(1)],[fixY targY bcitarget(2)],[r r 0; r r 0; bcitolerance]',[winColors;winColors;bciColors]);
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
        
       % msg('set 4 oval 0 %i %i %i %i %i %i',[0 0 cursorsize 0 255 0]);
    %    msg('set 5 oval 0 %i %i %i %i %i %i',[0 0 cursorsize 0 255 0]);
       
    %    msg('obj_switch 4 -5')
        msgAndWait(['diode 4']);
        msgAndWait(['diode 4']);
        %msgAndWait(['obj_on 4']); 
        %msg('obj_on 4')
        objon = 4;
        sendCode(codes.STIM1_ON);
    end
    thisStart = tic;
    
    if nargin >= 4 
      
         bciSuccess = 0;
         successclockflag = 0;
         freqrewardflag = 0;
         freqvals = [0 0];
         timeontargetcounter = 0;
         checkbcisuccessflag = 0;
        bciupdatetime = GetSecs;
        cursorupdatetime = bciupdatetime + cursoroffsettime;
        bciind =0;
        bciloopintervaltemp = bciloopinterval;
        sigflag = 0;
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
            
           
            %%%%%%%%%%%%%%%%%%% BCI code   %%%%%%%%%%%%%%%%%%%%%
            freqevents = [];
            if GetSecs - bciupdatetime >= bciloopintervaltemp
                 if matlabUDP2('check',sockets(2))
                     sigflag = 0;
                    %bciupdatetime = GetSecs;
                    while matlabUDP2('check',sockets(2))
                        freqevents = [freqevents; str2double(strsplit(matlabUDP2('receive',sockets(2))))];
                        %timevec = [timevec; (toc(thisStart)*1000)];
                        freqvals = [freqvals; freqevents];
                        sendCode(codes.SOUND_CHANGE);
                    end 
                    bciCursorTraj = freqvals;
                 else
                     if sigflag == 0 || (GetSecs - sigtime)>0.002
                        sigtime= GetSecs;
                        sigflag = 1;
                     end
                     
                 end
                 if ~isempty(freqevents)
                      if cursorvisibleflag == 1
                                msgclock = tic;
                                sendCode(codes.STIM1_ON)
                                msgAndWait('dset 4 oval 0 %i %i %i %i %i %i',[freqevents(end,1) freqevents(end,2) cursorsize cursorcolor(1) cursorcolor(2) cursorcolor(3)]);
                                sendCode(codes.STIM2_ON)
                                msgtime = toc(msgclock);
                                checkbcisuccessflag = 1;
                                
                      else
                          msgclock = tic;
                          while GetSecs-bciupdatetime < (bciloopinterval + cursoroffsettime)
                          end
                          msgtime = toc(msgclock);
                      end
                      if msgtime > cursoroffsettime
                          bciloopintervaltemp = bciloopinterval - (msgtime - cursoroffsettime);
                      else
                          bciloopintervaltemp = bciloopinterval;
                      end
                      bciupdatetime = GetSecs;   
                      timevec = [timevec; (toc(thisStart)*1000)];
                      sentcursors = [sentcursors;freqevents(end,1) freqevents(end,2) toc(thisStart)]; 
                 end
            end
             if checkbcisuccessflag == 1 
                 %bcifun(freqvals(end,:),bcitarget)
                    if ~isempty(freqvals)  && bcifun(freqvals(end,:),bcitarget)==1
                        if bcirewardflag == 1 && timeontargetcounter >= updatesintargetwin
                            fixSuccess = 1;
                            bciSuccess = 1;
                            break
                        else
                            timeontargetcounter = timeontargetcounter+1;
                        end
                        checkbcisuccessflag = 0;
                    else
                        checkbcisuccessflag = 0;
                        timeontargetcounter = 0;
                    end
             end
            %%%%%%%%%%%%%%%%%%%%%% End BCI Code %%%%%%%%%%%%%%%%%%%%%%%%%%5
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
end
