function result = ex_DistanceBCIBlock2BCI(e)
% ex file: ex_SaccadeTask
%
% Uses codes in the 2000s range to indicate stimulus types
% 2001 - visually guided saccade
% 2002 - memory guided saccade
% 2003 - delayed visually guided saccade
%
% General file for memory and visually guided saccade tasks 
%
% XML REQUIREMENTS
% angle: angle of the target dot from the fixation point 0-360
% distance: the distance of the target in pixels from the fixation point
% size: the size of the target in pixels
% targetColor: a 3 element [R G B] vector for the target color
% timeToFix: time in ms for initial fixation
% noFixTimeout: timeout punishment for aborted trial (ms)
% targetOnsetDelay: time after fixation before target appears
% fixDuration: length of initial fixation required
% targetDuration: duration that target is on screen
% stayOnTarget: length of target fixation required
% saccadeInitiate: maximum time allowed to leave fixation window
% saccadeTime: maximum time allowed to reach target
%
% Modified:
%
% 2012/10/22 by Matt Smith - update to helperTargetColor and to make it
% work when multiple 'trials' are passed in 'e'.
%
% 2012/11/09 by Matt Smith - update to allow "extraBorder" as an optional
% parameter in the XML file
%
% 2014/12/10 by Matt Smith - trying to consolidate some stuff.
% InterTrialPause put in at the end
%
% 2015/08/14 by Matt Smith - added ACQUIRE_TARG code before stayOnTarget
% window
%
% 2015/08/14 by Matt Smith - added recentering when detecting saccade
%



% %I think this should work, -ACS24Apr2012 %added the  cue for this particular trial (cue from xml just sets first miniblock) -ACS07oct2015 -added isCueTrial acs12oct2015



% xml: e.miniblocksize e.nCueTrials
% angle

%% initialize globals
    global params codes behav sockets bciCursorTraj trialData allCodes;
    bciCursorTraj = [];

%% Shake Hands with BCI Computer    
    handshakeflag = 0;
    
    if e.errorBadHandshakeFlag    
        matlabUDP2('send',sockets(2),'handshake');
        starthand = tic;
        while toc(starthand)<1
            if matlabUDP2('check',sockets(2))
                temp = matlabUDP2('receive',sockets(2));
                if str2double(temp)==1
                    handshakeflag = 1; 
                    break
                end
            end
        end

        if handshakeflag == 0 
            %error('BCI Computer Not Responding')
            error('exFunction:bci_aborted','BCI computer gave a bad handshake flag');
            % Don't need to send the code here because runex does it
            %sendCode(codes.BCI_ABORT);
        end
    end
    %prestim = tic;
    
    matlabUDP2('send',sockets(2),'trialend0');% ADDED 2/1/2019 RW: send this to ensure that params from previous trial are cleared. 

    
%% Initialize behav file    (FIX)
    if ~isfield(behav,'bcicorrect')
        thisdate = datevec(e(1).startBCIdate);
        behav.startBCI = mod(thisdate(3),2)+1;
        % automatically change start bci ind everyday. 
        behav.currBCI = behav.startBCI;
        % Later code will flip currBCI in first block, so need to update
        % actual start bci to match what will be used in first bci block
        behav.startBCI = (e(1).numbcis-behav.startBCI)+1; 
        behav.cursorloc = [];
        behav.bcicorrect = -1;
        behav.trialNum = 0;
        behav.cue = e(1).angle(randi(length(e(1).angle),1));
       % behav.superBlockCues = behav.cue;        
        behav.superblockcount = 0;
        behav.shamsleft = [];
        behav.bcisleft = [];
        %behav.flipFlag = 0;
        if isfield(e(1),'iterativeRecal') && e(1).iterativeRecal==1
            behav.superblockid = 0;
            behav.currRecal = 1;
        else
            behav.superblockid = 1;
        end
    else
        behav(end+1).bcicorrect = -1;
        behav(end).startBCI = behav(end-1).startBCI;
        behav(end).currBCI = behav(end-1).currBCI;
        behav(end).trialNum = behav(end-1).trialNum;
       % behav(end).superBlockCues = behav(end-1).superBlockCues;
       % behav(end).flipFlag = behav(end-1).flipFlag;
       % behav(end).cue = behav(end-1).cue;
        behav(end).superblockid = behav(end-1).superblockid;
        behav(end).superblockcount = behav(end-1).superblockcount;
        behav(end).shamsleft = behav(end-1).shamsleft;
        behav(end).bcisleft = behav(end-1).bcisleft;
        if isfield(e,'iterativeRecal') && e.iterativeRecal==1
            behav(end).currRecal = behav(end-1).currRecal;
        end
    end
%     if ~isempty(behav(1).bcicorrect)
%         trialData{6} = ['BCI Correct', num2str(sum([behav.bcicorrect]==1))];
%         trialData{7} = ['BCI Missed', num2str(sum([behav.bcicorrect]==0))];
%     end

    e = e(1); %in case more than one 'trial' is passed at a time...
%     recalflag = 0;
%     if isfield(e,'iterativeRecal') && e.iterativeRecal == 1 && behav(end).currRecal<=length(e.recalTrial)
%         recalflag = 1;
%     end

%% Determine block Identity (FIX)
% new code

% step 1: determine super block ID
% variables:
% e.bciinitblocklength
% e.bciblocklength
% e.shamblocklength
% behav(end).superblockid
% behav(end).superblockcount
switch behav(end).superblockid
    case 0       
        curblocklength = e.recalTrial(behav(end).currRecal);
        nextblock = 1;
    case 1; curblocklength = e.bciinitblocklength;nextblock = 2;% bciinitblock       
    case 2; curblocklength = e.bciblocklength;nextblock = 3;% bciblock
    case 3; curblocklength = e.shamblocklength;nextblock = 2;% shamblock
end
newblockflag = 0;
if behav(end).superblockcount >= curblocklength
    if behav(end).superblockid == 0
        behav(end).currRecal = behav(end).currRecal + 1;
        if behav(end).currRecal > length(e.recalTrial)
            behav(end).superblockid = nextblock;
            behav(end).superblockcount = 0;
            newblockflag = 1;
        else
            behav(end).superblockcount = 0;
        end
    else
        behav(end).superblockid = nextblock;
        behav(end).superblockcount = 0;
        newblockflag = 1;
    end
end
% step 2: determine miniblock values
% e.bciminiblocktrials; [1 1 1 1 1 1 1 1 1 0]
% e.shamminiblocktrials: [0 0 0 0 0 0 0 0 0 0]
% e.shampossible; [1 2 3 4 5]
% behav(end).shamsleft
% behav(end).bcisleft
if newblockflag == 1 || isempty(behav(end).bcisleft)% reset bcisleft and shamsleft in new block (case when superblocklength is not a multiple of miniblock length)
    switch behav(end).superblockid % if only one condition, do not need vector, ok to reset on every trial
        case 0; behav(end).bcisleft = 0;
        case 1
            behav(end).bcisleft = 1;
            if newblockflag == 1
                behav(end).currBCI = (e.numbcis-behav(end).currBCI)+1;% bciinitblock; set to bci 1 for now   
            end   
        case 2; behav(end).bcisleft = e.bciminiblocktrials;% bciblock
        case 3
            behav(end).bcisleft = e.shamminiblocktrials;
            if newblockflag == 1
                behav(end).currBCI = (e.numbcis-behav(end).currBCI)+1;% shamblock, may reset every time
            end
    end
end
if newblockflag == 1 || isempty(behav(end).shamsleft)
    behav(end).shamsleft = e.shampossible;
end

currtrialind = randi(length(behav(end).bcisleft));
currtrialid = behav(end).bcisleft(currtrialind);
newbcisleft = behav(end).bcisleft;
newbcisleft(currtrialind) = [];
bcinum = behav(end).currBCI;
% step 3: determine current trial ID
recalflag = 0;
if currtrialid == 0 && behav(end).superblockid~=0
    shamnumind = randi(length(behav(end).shamsleft));
    shamnum = behav(end).shamsleft(shamnumind);
    newshamsleft = behav(end).shamsleft;
    newshamsleft(shamnumind) = [];
    %bcinum = behav.currBCI; % always need this to be a valid bci index
    matlabUDP2('send',sockets(2),['sham',num2str(shamnum)]);
    % send message identifying sham trial num
elseif behav(end).superblockid==0
    recalflag = 1;
    shamnum = -1;
    %bcinum = behav.currBCI;; % set this to one so bci doesn't crash (will use cross talk modelparams file)
else
    shamnum = -1;
    %bcinum = behav.currBCI;
    matlabUDP2('send',sockets(2),['bcispec',num2str(bcinum)]);
    % send message identifying bci number
end
sendCode(1001+shamnum)
sendStruct(struct('bciTrial',e.bciTrial,'trialNum',behav(end).trialNum,'shamnum',shamnum,'bcinum',bcinum,'superblockID',behav(end).superblockid ));

% to do:
% modify bci code to hand multiple bci param files
% DONE: modify bci code to take in sham num and bci num
% DONE-ISH: clean up task code so old variables are handled correctly
% update xml file with new variables and new filecall
% make .mat file for bci playback
% DONE: add iterative recal block back
% test code

% old code
    objID = 2;
%     switch mod(behav(end).trialNum,e.miniblocksize),
%         case 0, %indicates start of new miniblock
%             if e.miniblocksize == 1
%                 behav(end).flipFlag = 0;
%             end
%             if behav(end).trialNum>0&&~behav(end).flipFlag, %don't flip before the first miniblock
%                 if length(behav(end).superBlockCues)==length(e.angle)
%                     behav(end).superBlockCues = [];
%                 end
%                 diff = setdiff(e.angle,behav(end).superBlockCues); %set the current cue to the "other" one
%                 thisPerm = randperm(length(diff),1);
%                 behav(end).cue = diff(thisPerm);
%                 behav(end).superBlockCues = [behav(end).superBlockCues behav(end).cue];
%                 behav(end).flipFlag = 1; %mark that flip has been done (in case the subject breaks fixation, etc. and we have to try the first trial again)
%             end;
%         case 1, %first trial after flip
%             behav(end).flipFlag = 0; %signal that flip will be needed next miniblock
%     end;
%     
%     cueTrialFlag = 0;
    targOnFlag = 0;
%     if mod(behav(end).trialNum,e.miniblocksize)<e.nCueTrials
%         if e.cueBCI==0 
%             e.bciTrial = 0;
%         end
%         targOnFlag = 1;
%         cueTrialFlag = 1;
%     end;
%     e.angle = (behav(end).cue); %pull the "real" cue from the behav structure
%     
    if e.convergeAnnulus == 1 && e.bciTrial == 1
        e.bciTrial = 0;
        stopbciflag = 1;
    else
        stopbciflag = 0;
    end    

%% Check other task input flags    
    result = 0;
    if e.showDelayTargetFlag == 1
        targetColorOther = e.targetColorOther;
    else
        targetColorOther = e.bgColor;
    end
    
    % flags for partial bci visibility
    if e.bciTrial == 0 
        cursorvisibleflag = 0;
        bcirewardflag = 0;
    else
        cursorvisibleflag = e.cursorVisible;
        bcirewardflag = e.bciRewardEnable;
    end
    
    if e.adeptImageFlag==1
        if e.adeptCloseorFar == 1
                        %e.angle
            %find(e.delayAngle==e.angle,1)
            thisindex =(find(e.delayAngle==e.angle,1));
            %thisindex
            newSuffix = e.adeptsuffix(e.indextouse(thisindex));
            %newSuffix
        else

            %length(e.indextouse)
            thisindex = e.indextouse(find(e.delayAngle==e.angle,1)+length(e.delayAngle));
            
            newSuffix = e.adeptsuffix(thisindex);
        end
    end
    
%% define bci function
    if isfield(e,'distfunflag')&&e.distfunflag == 1
        bcifun = @(x,y)(x<= y);
    else
    bcifun = @(x,y)(find(x==max(x))==y);
    end
   
    
    % take radius and angle and figure out x/y for saccade direction
    theta = deg2rad(e.angle);
    newX = round(e.distance*cos(theta));
    newY = round(e.distance*sin(theta));
%     while matlabUDP2('check',sockets(2))
%         dump = matlabUDP2('receive',sockets(2));
%     end 
%     matlabUDP2('send',sockets(2),['angle', num2str(e.angle)]);
%     matlabUDP2('send',sockets(2),['distance', num2str(e.distance)]);
%    matlabUDP2('send',sockets(2),['ncuetrials', num2str(e.nCueTrials)]);
    % now figure out if you need to shift the fixation point around so the
    % saccade will fit on the screen (e.g., for an 'amp' series). The
    % "extraborder" keeps the dot from ever getting within that many pixels
    % of the edge of the screen
    
%% some relic features of dirmem
    if isfield(e,'extraBorder') 
        extraborder = e.extraBorder; % use XML file if it's there
    else
        extraborder = 10; % default to 10 pixels
    end
    
    if (abs(newX) + e.size > (params.displayWidth/2 - extraborder))
        %disp('X exceeds limit, moving fix pt');
        shiftX = abs(newX) + e.size - params.displayWidth/2 + extraborder;
        if newX > 0
            e.fixX = e.fixX - shiftX;
            newX = newX - shiftX;
        else
            e.fixX = e.fixX + shiftX;
            newX = newX + shiftX;
        end
    end
    if (abs(newY) + e.size > (params.displayHeight/2 - extraborder))
        %disp('Y exceeds limit, moving fix pt');
        shiftY = abs(newY) + e.size - params.displayHeight/2 + extraborder;
        if newY > 0
            e.fixY = e.fixY - shiftY;
            newY = newY - shiftY;
        else
            e.fixY = e.fixY + shiftY;
            newY = newY + shiftY;
        end
    end
    
    % non targets positions. To do: add shift if exceed limits
    nonTargs = setdiff(e.allAngles,e.angle);
    nonTargX = zeros(length(nonTargs),1);
    nonTargY = zeros(length(nonTargs),1);
    for n = 1:length(nonTargs)
        localtheta = deg2rad(nonTargs(n));
        nonTargX(n) = round(e.distance*cos(localtheta));
        nonTargY(n) = round(e.distance*sin(localtheta));
    end
    
    
    
    % obj 1 is fix pt, obj 2 is target, diode attached to obj 2
    
%% set up the stimuli
    msg('set 1 oval 0 %i %i %i %i %i %i',[e.fixX e.fixY e.fixRad e.fixColor(1) e.fixColor(2) e.fixColor(3)]);
    if e.gratingsordots == 0
       % test = 1
        msg('set 2 oval 0 %i %i %i %i %i %i',[newX newY e.size e.targetColor(1) e.targetColor(2) e.targetColor(3)]);
    elseif e.gratingsordots == 1
       % test = 2
        msg('set 2 grating %i %f %f %f %f %i %i %i %f %f',[0  e.gratori(find(e.delayAngle==e.angle,1)) 0 e.gratspatial 0 newX newY e.size e.gratcontrast]);
    elseif e.gratingsordots == 2
       % test = 3
        runLine = e(1).runline;
        runString = '';
        while ~isempty(runLine)
            [tok, runLine] = strtok(runLine);

            while ~isempty(tok)
                [thisTok, tok] = strtok(tok,',');
                if strcmp('centerx',thisTok)
                    runString = [runString num2str(newX)];
                elseif strcmp('centery',thisTok)
                    runString = [runString num2str(newY)];
                elseif e.adeptImageFlag ==1 && strcmp('suffix',thisTok)
                    runString = [runString num2str(newSuffix)];
                else
                    runString = [runString num2str(eval(['e(1).' thisTok]))];
                end
            end

            runString = [runString ' '];
        end
        runString = [e(1).type ' ' runString(1:end-1)];
%         disp(['set ' num2str(objID) ' ' runString]);

        msg(['set ' num2str(2) ' ' runString]);
    end
    if isfield(e,'helperTargetColor')
        if isfield(e,'delayAngle') && isfield(e,'helperTargetColorAllOn')
            msg('set 3 oval 0 %i %i %i %i %i %i',[newX newY e.size e.helperTargetColorAllOn(1) e.helperTargetColorAllOn(2) e.helperTargetColorAllOn(3)]);
        else
            msg('set 3 oval 0 %i %i %i %i %i %i',[newX newY e.size e.helperTargetColor(1) e.helperTargetColor(2) e.helperTargetColor(3)]);
        end
    end
   % msg('set 4 oval 0 %i %i %i %i %i %i',[0 0 e.cursorSize e.cursorColor(1) e.cursorColor(2) e.cursorColor(3)]);
    %targetPost = 1-1/length(e.delayAngle);
%% set up feedback stim (annulus)    
    outerrad = round(e.minAnnulusSize+e.startAnnulusRad*(e.maxAnnulusSize-e.minAnnulusSize));
    innerrad = outerrad - e.annulusThickness;
    if  e.bciTrial == 1 || e.annulusvisibleflag == 1
        msg('set %i oval 0 %i %i %i %i %i %i',[4 e.fixX e.fixY (innerrad) e.bgColor(1) e.bgColor(2) e.bgColor(3)]);
        msg('set %i oval 0 %i %i %i %i %i %i',[5 e.fixX e.fixY (outerrad) e.cursorColor(1) e.cursorColor(2) e.cursorColor(3)]);
        
    end
%% set up multiple targets on the screen, if flag is set
    minnewobj = 6;
    numcolors = 100;
    if isfield(e,'delayAngle') && ~isempty(e.delayAngle)
        delaypositions = zeros(length(e.delayAngle),2);
        for stimind = 1:length(e.delayAngle)
                localtheta = deg2rad(e.delayAngle(stimind));
                localnewX = round(e.distance*cos(localtheta));
                delaypositions(stimind,1) = localnewX;
                localnewY = round(e.distance*sin(localtheta));
                delaypositions(stimind,2) = localnewY;
                if e.gratingsordotsother == 0
                    msg(['set ',num2str(stimind+minnewobj),' oval 0 %i %i %i %i %i %i'],[localnewX localnewY e.sizeOther targetColorOther(1) targetColorOther(2) targetColorOther(3)]);
                elseif e.gratingsordotsother == 1
                    msg(['set ',num2str(stimind+minnewobj),' grating %i %f %f %f %f %i %i %i %f %f'],[0  e.gratori(stimind) 0 e.gratspatial 0 localnewX localnewY e.sizeOther e.gratcontrast]);        
                else
                    runLine = e(1).runline;
                    runString = '';
                    while ~isempty(runLine)
                        [tok, runLine] = strtok(runLine);

                        while ~isempty(tok)
                            [thisTok, tok] = strtok(tok,',');
                            if strcmp('centerx',thisTok)
                                runString = [runString num2str(localnewX)];
                            elseif strcmp('centery',thisTok)
                                runString = [runString num2str(localnewY)];
                            else
                                runString = [runString num2str(eval(['e(1).' thisTok]))];
                            end
                        end

                        runString = [runString ' '];
                    end
                    runString = [e(1).type ' ' runString(1:end-1)];
            %         disp(['set ' num2str(objID) ' ' runString]);

                    msg(['set ' num2str(stimind+minnewobj) ' ' runString]);
                end
       end
    end
    bciTargRad = [params.targWinRad];
%% define the location of feedback (this may not work anymore...)
    switch e.feedbackType
        case 'target'
            cursorrad = [e.sizeOther];
            rgbvec = interpolateRGB2Gray(e.cursorColor,targetColorOther(1),linspace(0,1,numcolors),0);
        case 'fixation'
            cursorrad = [e.fixRad];
            rgbvec = interpolateRGB2Gray(e.cursorColor,e.fixColor,linspace(0,1,numcolors),0);
        otherwise
            cursorrad = [e.sizeOther];
            rgbvec = interpolateRGB2Gray(e.cursorColor,targetColorOther(1),linspace(0,1,numcolors),0);
    end
    msg(['diode ' num2str(objID)]);    
    
%     msgAndWait('ack'); %commented out 03Apr2013, seemed to be causing problems....
%     if cursorvisibleflag == 1
%         msg('obj_on 4'); % turn bci cursor on at the same time as fixation dot
%     end

%% start the task
    msgAndWait('obj_on 1');
    sendCode(codes.FIX_ON);
    
    if ~waitForFixation(e.timeToFix,e.fixX,e.fixY,params.fixWinRad);
        % failed to achieve fixation
        matlabUDP2('send',sockets(2),'trialend0');
        sendCode(codes.IGNORED);
        msgAndWait('all_off');
        sendCode(codes.FIX_OFF);
        waitForMS(e.noFixTimeout);
        result = codes.IGNORED;
        return;
    end
    sendCode(codes.FIXATE);
    if isfield(e,'fixJuice')
        if rand < e.fixJuice, giveJuice(1); end;
    end
   % predelaytime = toc(prestim)
    if ~waitForMS(e.targetOnsetDelay,e.fixX,e.fixY,params.fixWinRad)
        % hold fixation before stimulus comes on
        matlabUDP2('send',sockets(2),'trialend0');
        sendCode(codes.BROKE_FIX);
        msgAndWait('all_off');
        sendCode(codes.FIX_OFF);
        waitForMS(e.noFixTimeout);
        result = codes.BROKE_FIX;
        return;
    end
    if isfield(e,'delayAngle') && ~isempty(e.delayAngle)
        if strcmp(e.feedbackType,'annulus')
            msgAndWait(['obj_switch ',num2str((1:length(e.delayAngle))+minnewobj),' 4 5'])
        else
            msgAndWait(['obj_switch ',num2str((1:length(e.delayAngle))+minnewobj),' -4 -5'])
        end
    end
    % Decision point - is this VisGuided, Delay-VisGuided, or Mem-Guided
    if (e.targetOnsetDelay == e.fixDuration)
        % Visually Guided Saccade
        sendCode(2001); % send code specific to this stimulus type
        % turn fix pt off and target on simultaneously
        msgAndWait('obj_switch 2 -1');
%         msg('queue_begin');
%         msg('obj_on 2');
%         msg('obj_off 1'); 
%         msgAndWait('queue_end');
        sendCode(codes.FIX_OFF);
        sendCode(codes.TARG_ON);
    elseif ((e.targetOnsetDelay + e.targetDuration) < e.fixDuration) 
        % Memory Guided Saccade
        sendCode(2002); % send code specific to this stimulus type
        if targOnFlag == 1 
            msgAndWait('obj_on 2');
            sendCode(codes.TARG_ON);
        end
        matlabUDP2('send',sockets(2),'pretrialstart');
        if ~waitForMS(e.targetDuration,e.fixX,e.fixY,params.fixWinRad)
            % didn't hold fixation during target display
            matlabUDP2('send',sockets(2),'trialend0');
            sendCode(codes.BROKE_FIX);
            msgAndWait('all_off');
            sendCode(codes.TARG_OFF);
            sendCode(codes.FIX_OFF);
            waitForMS(e.noFixTimeout);
            result = codes.BROKE_FIX;
            return;
        end     

            msgAndWait('obj_off 2');
            sendCode(codes.TARG_OFF);
            prewait = tic;
            msgAndWait('diode 4');
%         if ~waitForMS(e.postTargDelay,e.fixX,e.fixY,params.fixWinRad)
%             % didn't hold fixation during target display
%             sendCode(codes.BROKE_FIX);
%             msgAndWait('all_off');
%             sendCode(codes.FIX_OFF);
%             waitForMS(e.noFixTimeout);
%             result = codes.BROKE_FIX;
%             return;
%         end     
        %
        if ~waitForMS(e.preBCIDelay-toc(prewait)*1000,e.fixX,e.fixY,params.fixWinRad)
            % didn't hold fixation during prebci delay
            matlabUDP2('send',sockets(2),'trialend0');
            sendCode(codes.BROKE_FIX);
            msgAndWait('all_off');
            sendCode(codes.FIX_OFF);
            waitForMS(e.noFixTimeout);
            result = codes.BROKE_FIX;
            return;
        end    
        
        sendCode(codes.ALIGN);
        matlabUDP2('send',sockets(2),['trialstart', num2str(numel(allCodes))]);
        
        %testwait = tic;
        if ~waitForMS(e.postBCIDelay,e.fixX,e.fixY,params.fixWinRad)
            % didn't hold fixation during postbcidelay
            matlabUDP2('send',sockets(2),'trialend0');
            sendCode(codes.BROKE_FIX);
            msgAndWait('all_off');
            sendCode(codes.FIX_OFF);
            waitForMS(e.noFixTimeout);
            result = codes.BROKE_FIX;
            return;
        end    
        waitRemainder = e.fixDuration - (e.targetOnsetDelay + e.targetDuration+e.postBCIDelay + e.preBCIDelay);
        [fixsuccess, bcisuccess, behav(end).cursorloc,timevec] = waitForMSmatlabUDP_distanceBCI(waitRemainder,e.fixX,e.fixY,params.fixWinRad,newX,newY,find(e.angle == e.delayAngle),bciTargRad,bcifun,...
                                                                        e.updatesOnTarget,cursorvisibleflag,rgbvec,numcolors,cursorrad,bcirewardflag,minnewobj+(1:length(e.delayAngle)),delaypositions,e.feedbackType,e);
        if recalflag == 1
            if isfield(e,'recalMissReward') && e.recalMissReward == 1 && fixsuccess==1
                bciSuccess = 1;
            end
        end
                                                                %         waittime = toc(testwait);
%         %display([newX newY])
%         %display(behav(end).cursorloc(end,:))
%        numfreqs = size(behav(end).cursorloc,1);
 %        display(behav(end).cursorloc)
 %        display(diff([0;timevec]))
%         if (waitRemainder/50)-1>numfreqs
%             display(waittime)
%             display(waitRemainder)
%             display(numfreqs)
%             display(behav(end).cursorloc)
             %display(diff([0;timevec]))
%         end
        
        if fixsuccess ~= 1
            % didn't hold fixation during period after target offset
            matlabUDP2('send',sockets(2),'trialend0');
            sendCode(codes.BROKE_FIX);
            msgAndWait('all_off');
            sendCode(codes.FIX_OFF);
            waitForMS(e.noFixTimeoutBCI);
            result = codes.BROKE_FIX;            
            return;
        elseif bcisuccess==1

            if e.bciTrial == 1
                if bcirewardflag == 1
                    matlabUDP2('send',sockets(2),'trialend2');
                    sendCode(codes.BCI_CORRECT);
                    result = codes.BCI_CORRECT;
                    %sendCode(codes.CORRECT);
                    if isfield(e,'delayAngle') && ~isempty(e.delayAngle)
                        msgAndWait(['obj_switch ',num2str(-1*((1:length(e.delayAngle))+minnewobj)), ' -4 -5'])
                    end
                    sendCode(codes.REWARD);
                    giveJuice(e.bciJuice);
                    waitForMS(e.BCIPause);
                    msgAndWait('all_off');
                    sendCode(codes.FIX_OFF);
                    behav(end).bcicorrect = 1;
                    behav(end).trialNum = behav(end).trialNum + 1;
                    if shamnum ~= -1
                        behav(end).shamsleft = newshamsleft;
                    end
                    behav(end).bcisleft = newbcisleft;
                    behav(end).superblockcount = behav(end).superblockcount + 1;
                    if recalflag == 1
                        totalbcitrials = sum([behav.bcicorrect]==1) + sum([behav.bcicorrect]==0);
                        if (totalbcitrials==e.recalTrial(behav(end).currRecal))
                            matlabUDP2('send',sockets(2),'recalibrate');
                            %behav(end).currRecal = behav(end).currRecal + 1;
                            waitForMS(e.calibrationpause);
                        end
                    end
                    return
                else
                    sendCode(codes.BCI_CORRECT);
                    behav(end).bcicorrect =  1;                
                end
            else
                if stopbciflag == 1
                    matlabUDP2('send',sockets(2),'trialend0');
                end
                behav(end).bcicorrect = -1;
            end

        elseif bcisuccess == 0
            if e.bciTrial == 1
                matlabUDP2('send',sockets(2),'trialend1');
                sendCode(codes.BCI_MISSED);
                result = codes.BCI_MISSED;
                behav(end).bcicorrect = 0;
                if shamnum ~= -1
                    behav(end).shamsleft = newshamsleft;
                end
                behav(end).bcisleft = newbcisleft;
                behav(end).superblockcount = behav(end).superblockcount + 1;
            else
                if stopbciflag == 1
                    matlabUDP2('send',sockets(2),'trialend0');
                end
                behav(end).bcicorrect = -1;
            end
        end
        %msgAndWait('obj_off 1'); 
        msgAndWait('obj_off 1 4 5'); 
        sendCode(codes.FIX_OFF);
        minrxn = tic;       
    else
        warning('*** EX_SACCADETASK: Condition not valid');
        %%% should there be some other behavior here?
        return;
    end
    
    msgAndWait('all_off');
%     sendCode(codes.CORRECT);
    sendCode(codes.TARG_OFF);
    sendCode(codes.REWARD);
    if e.bciTrial == 0 || e.bciRewardEnable==0
        result = codes.CORRECT;
        sendCode(codes.CORRECT);
    end
    if  e.bciTrial == 0 || e.bciRewardEnable==0
        giveJuice(e.calibJuice);
        behav(end).trialNum = behav(end).trialNum + 1;
    end
    if recalflag == 1
        totalbcitrials = sum([behav.bcicorrect]==1) + sum([behav.bcicorrect]==0);
        if (totalbcitrials==e.recalTrial(behav(end).currRecal))
            matlabUDP2('send',sockets(2),'recalibrate');
            %behav(end).currRecal = behav(end).currRecal + 1;
            waitForMS(e.calibrationpause);
        end
    end

    if isfield(e,'InterTrialPause')
        waitForMS(e.InterTrialPause); %this was for Wile E to lengthen time between trials SBK
    end
    
