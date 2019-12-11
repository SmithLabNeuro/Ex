function [fixSuccess, bciSuccess, sentcursors,timevec ]= waitForMSmatlabUDP_distanceBCI(waitTime,fixX,fixY,r,targX,targY,bcitarget, bcitolerance,bcifun,updatesintargetwin,cursorvisibleflag,cursorcolor,numcolors,cursorsize,bcirewardflag,ininds,delaytargpos,feedbackType,e,varargin)
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
% 2016/5/24 Ryan Williamson added functionality for communicating with a
% bci system

global params;
global codes;
global sockets;
global behav;
%global bciCursorTraj;
timevec = [];
% check inputs
complexAngles = exp(1i*deg2rad(e.delayAngle));
winColors = [255 255 0];
bciloopintervalpre = 0.04996; % this was empirically set to sync with the monitor updates which are not exactly every 10 ms
cursoroffsettime = 0.020;
distractorflag = 0;
bciloopinterval = bciloopintervalpre - cursoroffsettime;
recenterFlag = false;
sentcursors = [];
comparecolorvec = 1-linspace(0,1,numcolors); % create mapping from normalized LL to color with grey equal to 0
lowerbound = e.minAnnulusSize+e.lowerboundprop*(e.maxAnnulusSize-e.minAnnulusSize);
startrad = e.minAnnulusSize+e.startAnnulusRad*(e.maxAnnulusSize-e.minAnnulusSize);

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
    bciColors = [0 255 0];
    drawFixationWindows([fixX targX],[fixY targY],[r bcitolerance],[winColors; bciColors]);
    
elseif nargin ~=1
    error('waitForMS can have exactly 1, 4 or 5 input arguments');
end

fixSuccess = 1;

if cursorvisibleflag == 1
    sendCode(codes.STIM1_ON);
end
thisStart = tic;

if nargin >= 4
    
    bciSuccess = 0;
    freqvals = zeros(1,length(e.delayAngle));
    timeontargetcounter = 0;
    checkbcisuccessflag = 0;
    bciupdatetime = GetSecs;
    bciloopintervaltemp = bciloopinterval;
    sigflag = 0;
    while (toc(thisStart)*1000) <= waitTime
        loopTop = GetSecs;
        d=samp;
        
        eyePos = projectCalibration(d(end,:)); %changed to new project calibration function (supports polynomial regressors) -ACS 29Oct2013
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
        %distractorOnTime
        %distractorOn
        %distractorDuration
        %random: distractorOffset
        %distColor
        % define the threshold distance
        if length(e.mindist) == length(e.delayAngle)
            currmindist = e.mindist(bcitarget);
        elseif length(e.mindist) == 1
            currmindist = e.mindist;
        else
            error('e.mindist should be either length 1 or length equal to e.angle')
        end
        if isfield(e,'distractorOn')&&e.distractorOn == 1
            if (toc(thisStart)*1000)+e.preBCIDelay>=e.distractorOnTime&&distractorflag == 0
                msgAndWait('obj_on 6');
                sendCode(codes.TARG2_ON);
                distractorflag = 1;
            elseif (toc(thisStart)*1000)+e.preBCIDelay>=(e.distractorOnTime+e.distractorDuration)&&distractorflag == 1
                msgAndWait('obj_off 6');
                sendCode(codes.TARG2_OFF);
                distractorflag = 2;
            end
        end
        
        
        freqevents = [];
        if GetSecs - bciupdatetime >= bciloopintervaltemp
            if matlabUDP2('check',sockets(2))
                sigflag = 0;
                while matlabUDP2('check',sockets(2))
                    a = str2double(strsplit(matlabUDP2('receive',sockets(2))));
                    freqevents = [freqevents; a];
                    freqvals = [freqvals; a];
                    sendCode(codes.SOUND_CHANGE);
                end
            else
                if sigflag == 0 || (GetSecs - sigtime)>0.002
                    sigtime= GetSecs;
                    sigflag = 1;
                end       
            end
            
            if e.convergeAnnulus == 0 &&~isempty(freqevents)
                if cursorvisibleflag == 1
                    msgclock = tic;
                    sendCode(codes.STIM1_ON)
                    posterior = freqevents(end);
                    targetcolors = zeros(length(posterior),3);
                    for n = 1:size(targetcolors,1)
                        [~,I] = min(abs(comparecolorvec-posterior(n)));
                        targetcolors(n,:) = cursorcolor(I,:);
                    end
                    msgstring = [];
                    targetdist = posterior;
                    display('I am here but the lights are off')
                    targetdist = max(0,e.annulusSizeAtThresh*posterior);
%                     if isfield(e,'iterativeRecal') && e.iterativeRecal == 1 && behav(end).currRecal<=length(e.recalTrial)
%                         convergeTime = (waitTime-e.convergeDelay);
%                         if (toc(thisStart)*1000) <= convergeTime% freeze annulus T ms before end of delay
%                             delayRemaining = 1-(toc(thisStart)*1000)/convergeTime;
%                             %delayRemaining = 1-(toc(thisStart)*1000)/convergeTime;
% %                             if e.startAnnulusRad > e.lowerboundprop
% %                                 outerradnobci = round(e.minAnnulusSize+delayRemaining*(e.maxAnnulusSize-e.minAnnulusSize));
% %                             end
%                         end
%                         % 2018/11/28 RW make the annulus stay just above the
%                         % threshold so that he does not get out of converge
%                         % trials early
%                         targetdist = (targetdist*e.recalPercent(behav(end).currRecal) + max(e.annulusSizeAtThresh+0.0001,delayRemaining)*(1-e.recalPercent(behav(end).currRecal)));
%                         if isfield(e,'recalMissReward') && e.recalMissReward == 1
%                             bciSuccess = 1;
%                         end
%                     end
%                     
                     outerrad = round(e.minAnnulusSize+targetdist*(e.maxAnnulusSize-e.minAnnulusSize));
                    %display(targetdist)
                    %%%%%%% code for gradually transfering control to animal %%%%%

                  %%%%% end code for gradually transfering control to animal %%%%%
                    
                    
                    
                    
                    
                    innerrad = outerrad - e.annulusThickness;
                    msgstring1 = sprintf('mset %i oval 0 %i %i %i %i %i %i',[5 fixX fixY (outerrad) e.cursorColor(1) e.cursorColor(2) e.cursorColor(3)]);
                    msgstring2 = sprintf(' mset %i oval 0 %i %i %i %i %i %i',[4 fixX fixY (innerrad) e.bgColor(1) e.bgColor(2) e.bgColor(3)]);
                    msgstring = [msgstring1 msgstring2];
                    msgAndWait(msgstring);
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
                sentcursors = [sentcursors;targetdist toc(thisStart)];
            elseif e.convergeAnnulus == 1 && e.annulusvisibleflag==1
                convergeTime = (waitTime-e.convergeDelay);
                
                if (toc(thisStart)*1000) <= convergeTime% freeze annulus T ms before end of delay
                    if e.noisyAnnulusFlag==1&&(e.noisyFlatTime+e.preBCIDelay+e.postBCIDelay)<e.fixDuration
                        flattimeleft = (toc(thisStart)*1000);
                        delayRemaining = 1-(toc(thisStart)*1000-e.noisyFlatTime)/(convergeTime-e.noisyFlatTime);
                    else
                        delayRemaining = 1-(toc(thisStart)*1000)/convergeTime;
                    end
                    if e.noisyAnnulusFlag==1 && (e.noisyFlatTime)<convergeTime && flattimeleft < e.noisyFlatTime
                        outerrad=round(startrad+e.initSinusoidGain*sin(pi*(flattimeleft/(1000*bciloopintervalpre))/e.noisyFreq1)+e.initSinusoidGain*sin(pi*(flattimeleft/(1000*bciloopintervalpre))/e.noisyFreq2));
                        innerrad = outerrad - e.annulusThickness;
                    else
                        %delayRemaining = 1-(toc(thisStart)*1000)/convergeTime;
                        if e.startAnnulusRad>e.lowerboundprop
                            outerrad = round(lowerbound+delayRemaining*(startrad-lowerbound));
                            innerrad = outerrad - e.annulusThickness;
                        end
                    end
                end
                msgstring1 = sprintf('mset %i oval 0 %i %i %i %i %i %i',[5 fixX fixY (outerrad) e.cursorColor(1) e.cursorColor(2) e.cursorColor(3)]);
                msgstring2 = sprintf(' mset %i oval 0 %i %i %i %i %i %i',[4 fixX fixY (innerrad) e.bgColor(1) e.bgColor(2) e.bgColor(3)]);
                msgstring = [msgstring1 msgstring2];
                msgAndWait(msgstring);
                sendCode(codes.STIM2_ON)
            end
        end
        if checkbcisuccessflag == 1
            %bcifun(freqvals(end,:),bcitarget)
            checkbcisuccessflag = 0;

            bcifunout = bcifun(targetdist,e.annulusSizeAtThresh);% the posterior should always be normalized so that correct = 1
            if ~isempty(freqvals)  && length(bcifunout)==1 && bcifunout==1
                timeontargetcounter = timeontargetcounter+1;
                if bcirewardflag == 1 && timeontargetcounter >= updatesintargetwin
                    fixSuccess = 1;
                    bciSuccess = 1;
                    break
                end
            else
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
