function result = ex_driftchoice(e)
% ex file: ex_driftchoice
%
%
%
% Modified:
%
global params codes behav allCodes;

e = e(1); %in case more than one 'trial' is passed at a time...

backdoor = e(1).tooSlow;

%initialize behavior-related stuff:
if ~isfield(behav,'score')||~isfield(behav,'targAmp')||~isfield(behav,'targPickHistory')
    behav.score = [];
    behav.targAmp = sort(e.startTargAmp,'ascend');
    behav.targPickHistory = [];
    behav.trialNum = [];
    behav.RT = [];
    behav.showHelp = [];
end;

%Constant stuff:
frameMsec = params.displayFrameTime*1000;
stimulusDuration = e.stimulusDuration;
targetDuration = e.targetDuration;
waitAfterTarget = e.waitAfterTarget;
minShortTargetOnset = e.minShortTargetOnset;
maxShortTargetOnset = e.maxShortTargetOnset;
minLongTargetOnset = maxShortTargetOnset;
maxLongTargetOnset = stimulusDuration-targetDuration;
shortLongTrialProportion = e.shortLongTrialProportion;
cueDuration = e.cueDuration;
cueWin = @hann;
cueFs = 44100;
cueFreq = e.cueFreq;
startColors = [e.startColor1;e.startColor2];
endColors = [e.endColor1;e.endColor2];
orientations = [e.orientation1,e.orientation2];
isCatch = e.isCatch;
extraReward = e.extraReward;
if isfield(e,'helpFixProb'), showHelp = rand<e.helpFixProb; else showHelp = true; end;

%Variable stuff:
targetObject = abs(((1-e.isValid)*3)-e.cue);
isShortTrial = rand<=shortLongTrialProportion;
if isCatch
    targetOnset = stimulusDuration+1;
elseif isShortTrial
    targetOnset = randi(maxShortTargetOnset-minShortTargetOnset)+minShortTargetOnset;
else
    targetOnset = randi(maxLongTargetOnset-minLongTargetOnset)+minLongTargetOnset;
end;
targetPeak = targetOnset+round(targetDuration./2); %frames
targetEnd = targetOnset+targetDuration;

targAmpPick = e.targAmpPick;
targAmp = behav.targAmp(targAmpPick)*e.temporal;
sendStruct(struct('targAmp',targAmp)); %I think this should work, -ACS24Apr2012
sendStruct(struct('targOnset',targetOnset));

switch e.cueMap
    case 1 %spatial cue
        colorPick = [e.colorPick 3-e.colorPick];
        targColors = [startColors(colorPick(1),:);endColors(colorPick(1),:)];
        distColors = [startColors(colorPick(2),:);endColors(colorPick(2),:)];
        posPick = (1-targetObject)*2+1;
        targX = e.centerx*posPick;
        targY = e.centery;
        distX = -e.centerx*posPick;
        distY = e.centery;
        cueX = e.centerx*(-2*e.cue+3);
        cueY = e.centery;
        cueColor = fix([255 255 255]*e.visCueBrightness);
    case 2 %featural cue
        targColors = [startColors(targetObject,:);endColors(targetObject,:)];
        distColors = [startColors(3-targetObject,:);endColors(3-targetObject,:)];
        posPick = e.posPick;
        targX = e.centerx*posPick;
        targY = e.centery;
        distX = -e.centerx*posPick;
        distY = -e.centery;
        cueX = 0; cueY = 0;
        cueColor = fix(endColors(targetObject,:)*e.visCueBrightness);
    otherwise
        error('driftchoice:unknownCueMap','unknown cue mapping');
end;
orientations = orientations([e.oriPick 3-e.oriPick]); %randomize orientations

cueColor = fix([131 132 90]*e.visCueBrightness); %rough hard code for now -ACS 21May2014

msgAndWait('set 5 rgbgrating %i %f %f %f %f %i %i %i %f %i %i %i %i %i %i %i %f %i %i %i',...
    [stimulusDuration  orientations(1) e.phase e.spatial e.temporal targX targY e.radius e.contrast targColors(:)' e.cmapReso,...
    targAmp  targetOnset targetPeak targetEnd]);
msgAndWait('set 6 rgbgrating %i %f %f %f %f %i %i %i %f %i %i %i %i %i %i %i %f %i %i %i',...
    [stimulusDuration  orientations(2) e.phase e.spatial e.temporal distX distY e.radius e.contrast distColors(:)' e.cmapReso,...
    e.temporal targetOnset targetPeak targetEnd]);
msgAndWait('set 7 oval %i %i %i %i %i %i %i',[floor(cueDuration*1000/frameMsec) cueX cueY e.radius+e.margin 127 127 127]); %gray margin around RF
msgAndWait('set 8 rect %i %i %i %i %i %i %i %i',[floor(cueDuration*1000/frameMsec) -cueX cueY abs(cueX) 400 127 127 127]);
msgAndWait('set 9 oval %i %i %i %i %i %i %i',[floor(cueDuration*1000/frameMsec) cueX cueY e.visCueRad cueColor]);

msgAndWait('set 1 oval 0 %i %i %i %i %i %i',[e.fixX e.fixY e.fixRad 255 255 0]); %constant central fixation (yellow)
if e.showHoles
    if show,
        msgAndWait('set 2 oval 0 %i %i %i %i %i %i',[targX targY e.fixRad e.helpFixColor]); %target fixation (blue)
    else
        msgAndWait('set 2 oval 0 %i %i %i %i %i %i',[targX targY e.fixRad 127 127 127]); %target fixation (blue)
    end;
    msgAndWait('set 3 oval 0 %i %i %i %i %i %i',[targX targY e.fixRad 127 127 127]); %'hole' in target grating
    msgAndWait('set 4 oval 0 %i %i %i %i %i %i',[distX distY e.fixRad 127 127 127]); %'hole' in distracter grating
end;
msg('diode 5');

msgAndWait('ack');

msgAndWait('obj_on 1');
sendCode(codes.FIX_ON);

if ~waitForFixation(e.timeToFix,e.fixX,e.fixY,params.fixWinRad)
    % failed to achieve fixation
    sendCode(codes.IGNORED);
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    waitForMS(1000); % no full time-out in this case
    result = codes.IGNORED;
    return;
end

sendCode(codes.FIXATE);

if ~waitForMS(e.preStimFix,e.fixX,e.fixY,params.fixWinRad)
    % hold fixation before stimulus comes on
    sendCode(codes.BROKE_FIX);
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    waitForMS(e.noFixTimeout);
    result = codes.BROKE_FIX;
    return;
end

%the cue:
if cueDuration>0
    timebase = 0:1/cueFs:cueDuration;
    cueSound = feval(cueWin,numel(timebase))'.*sin(2.*pi.*cueFreq(e.cue).*timebase); %pure tone
%     cueSound = feval(cueWin,numel(timebase))'.*chirp(timebase,cueFreq(e.cue),cueDuration,cueFreq(3-e.cue)); %swept frequency chirp
    if e.cue==2, cueAmplitude = fliplr(e.cueAmplitude); else cueAmplitude = e.cueAmplitude; end;
    if isfield(params,'flipSpeakers')&&params.flipSpeakers==1, cueAmplitude = fliplr(cueAmplitude); end;
    cueSound = cueSound(:)*cueAmplitude;
    codeStr = sprintf('STIM%d_ON',e.cue);
    sound(cueSound,cueFs);
    if any(e.cueAmplitude>0)
        msgAndWait('queue_begin');
            msg('obj_on 7');
            msg('obj_on 8');
            msg('obj_on 9');
        msgAndWait('queue_end');
        sendCode(codes.(codeStr));
    else
        sendCode(codes.NO_STIM);
    end;
end;

if ~waitForMS(e.cueTargetInterval,e.fixX,e.fixY,params.fixWinRad)
    % failed to keep fixation
    sendCode(codes.BROKE_FIX);
    msgAndWait('all_off');
    sendCode(codes.STIM_OFF);
    sendCode(codes.FIX_OFF);
    waitForMS(e.noFixTimeout);
    result = codes.BROKE_FIX;
    return;
end

%intro
msgAndWait('queue_begin');
msg('obj_on 5');
msg('obj_on 6');
if e.showHoles
    msg('obj_on 4');
    msg('obj_on 3');
end;
msgAndWait('queue_end');
sendCode(codes.STIM_ON);

hitFlag = false;
winColors = [0,255,0;255,0,0];
%require hold fixation until target onset:
if ~waitForMS((targetOnset+waitAfterTarget).*frameMsec,e.fixX,e.fixY,params.fixWinRad)
    % failed to keep fixation
    choiceWin = waitForFixation(8.*frameMsec,[targX distX],[targY distY],params.targWinRad.*[1 1]); %note that target window is always '1'
    switch choiceWin
        case 1
            if e.isValid, choicePos = e.cue; else choicePos = 3-e.cue; end;
            result = codes.FALSEALARM;
        case 2
            if e.isValid, choicePos = 3-e.cue; else choicePos = e.cue; end;
            result = codes.FALSEALARM;
        otherwise
            choicePos = 0;
            result = codes.BROKE_FIX;
    end;
    sendCode(codes.(sprintf('CHOICE%d',choicePos)));
    msgAndWait('all_off');
    sendCode(codes.STIM_OFF);
    sendCode(codes.FIX_OFF);
    behav.score(end+1) = nan;
    behav.showHelp(end+1) = showHelp;
    behav.trialNum(end+1) = length(allCodes);
    behav.RT(end+1) = nan;
    adjustTarget(e,targAmpPick);
    waitForMS(e.noFixTimeout);
    return;
elseif isCatch %Correct withhold to a catch trial:
    msgAndWait('all_off');
    sendCode(codes.WITHHOLD);
    sendCode(codes.STIM_OFF);
    sendCode(codes.FIX_OFF);

    behav.score(end+1) = nan;
    behav.showHelp(end+1) = showHelp;
    behav.trialNum(end+1) = length(allCodes);
    behav.RT(end+1) = nan;
    adjustTarget(e,targAmpPick);
    giveJuice(2); % 2 for catch trial
    sendCode(codes.REWARD);
    result = codes.WITHHOLD;
    return;
else
    sendCode(codes.TARG_ON);
    targOnTime = tic;
end;

if e.showHoles
    choiceWin = waitForFixation((targetPeak-targetOnset).*frameMsec,[targX distX],[targY distY],params.targWinRad.*[1 1],winColors); %note that target window is always '1'
    switch choiceWin
        case 1 %saccade to target
            if ~waitForMS(e.stayOnTarget,targX,targY,params.targWinRad) %require to stay on for a while, in case the eye 'accidentally' travels through target window
                % failed to keep fixation
                sendCode(codes.BROKE_TARG);
                msgAndWait('all_off');
                sendCode(codes.STIM_OFF);
                sendCode(codes.FIX_OFF);
                behav.score(end+1) = nan;
                behav.showHelp(end+1) = showHelp;
                behav.trialNum(end+1) = length(allCodes);
                behav.RT(end+1) = nan;
                adjustTarget(e,targAmpPick);
                waitForMS(e.noFixTimeout);
                result = codes.BROKE_TARG;
                return;
            end;
            hitFlag = true;
            if toc(targOnTime)<=backdoor
                sendCode(codes.CORRECT);
            else
                sendCode(codes.MISSED);
            end;
            behav.score(end+1) = 1;
            behav.showHelp(end+1) = showHelp;
            behav.trialNum(end+1) = length(allCodes);
            behav.RT(end+1) = toc(targOnTime);
        case 2 %saccade to distracter
            if toc(targOnTime)<=backdoor
                sendCode(codes.WRONG_TARG);
                result = codes.WRONG_TARG;
            else
                sendCode(codes.MISSED);
                result = codes.MISSED;
            end;
            msgAndWait('all_off');
            sendCode(codes.STIM_OFF);
            sendCode(codes.FIX_OFF);
            behav.score(end+1) = 0;
            behav.showHelp(end+1) = showHelp;
            behav.trialNum(end+1) = length(allCodes);
            behav.RT(end+1) = nan;
            adjustTarget(e,targAmpPick);
            waitForMS(2.*e.noFixTimeout);
            return;
        otherwise
            %do nothing (yet)
    end;
    msg('obj_on 2');
    sendCode(codes.FIX_MOVE);
end

if ~hitFlag
    choiceWin = waitForFixation((stimulusDuration-targetOnset).*frameMsec,[targX distX],[targY distY],params.targWinRad.*[1 1],winColors(1:2,:)); %note that target window is always '1'
    %     fprintf('\nChoice:%i, Target:%i',choiceWin,targetObject); %debugging feedback
    switch choiceWin
        case 0
            % didn't reach target or distracter
            if waitForFixation(1,e.fixX,e.fixY,params.fixWinRad)
                sendCode(codes.NO_CHOICE);
                msgAndWait('all_off');
                sendCode(codes.STIM_OFF);
                sendCode(codes.FIX_OFF);
                behav.score(end+1) = 0;
                behav.showHelp(end+1) = showHelp;
                behav.trialNum(end+1) = length(allCodes);
                behav.RT(end+1) = nan;
                result = codes.NO_CHOICE;
            else
                msgAndWait('all_off');
                sendCode(codes.STIM_OFF);
                sendCode(codes.FIX_OFF);
                behav.score(end+1) = nan;
                behav.showHelp(end+1) = showHelp;
                behav.trialNum(end+1) = length(allCodes);
                behav.RT(end+1) = nan;
                result = codes.BROKE_FIX;
            end;
            adjustTarget(e,targAmpPick);
            return;
        case 1 %the target window
            if ~waitForMS(e.stayOnTarget,targX,targY,params.targWinRad) %require to stay on for a while, in case the eye 'accidentally' travels through target window
                % failed to keep fixation
                sendCode(codes.BROKE_TARG);
                msgAndWait('all_off');
                sendCode(codes.STIM_OFF);
                sendCode(codes.FIX_OFF);
                behav.score(end+1) = nan;
                behav.showHelp(end+1) = showHelp;
                behav.trialNum(end+1) = length(allCodes);
                behav.RT(end+1) = nan;
                adjustTarget(e,targAmpPick);
                waitForMS(e.noFixTimeout);
                result = codes.BROKE_TARG;
                return;
            end;
            if toc(targOnTime)<=backdoor
                sendCode(codes.CORRECT);
            else
                sendCode(codes.MISSED);
            end;                   
            behav.score(end+1) = 1;
            behav.showHelp(end+1) = showHelp;
            behav.trialNum(end+1) = length(allCodes);
            behav.RT(end+1) = toc(targOnTime);
        otherwise
            % incorrect choice
            %for a wrong choice, immediately score it, so that the subject
            %can't then switch to the other stimulus.
            if toc(targOnTime)<=backdoor
                sendCode(codes.WRONG_TARG); 
                result = codes.WRONG_TARG; %change the result code to indicate a wrong choice
            else
                sendCode(codes.MISSED);
                result = codes.MISSED;
            end;               
            msgAndWait('all_off');
            sendCode(codes.STIM_OFF);
            sendCode(codes.FIX_OFF);
            behav.score(end+1) = 0;
            behav.showHelp(end+1) = showHelp;
            behav.trialNum(end+1) = length(allCodes);
            behav.RT(end+1) = nan;
            adjustTarget(e,targAmpPick);
            waitForMS(2.*e.noFixTimeout);
            return;
    end;
end;

% turn off stimulus
msgAndWait('all_off');
sendCode(codes.STIM_OFF);
sendCode(codes.FIX_OFF);

adjustTarget(e,targAmpPick);
if toc(targOnTime)<=backdoor %new to close the back door on the response window 13May2013 -ACS
    if e.isValid==1
        sound(cueSound,cueFs);
        giveJuice(1+extraReward); % extra reward for valid corrects, 1 for invalid corrects
    else
        giveJuice(1);
    end
    sendCode(codes.REWARD);
    result = codes.CORRECT;
else
    behav.score(end) = 0;
    result = codes.MISSED;
end;

function adjustTarget(e,targAmpPick)
%Adjust target amplitude based on behavior (and update on-screen
%score):
global trialData behav
behav.targPickHistory(end+1) = targAmpPick;
if ~e.showHoles %don't adjust if you're just giving the subject the answer
    %         display(behav.targPickHistory);
    %         display(targAmpPick);
    if ~isnan(behav.score(end))&&sum(~isnan(behav.score(behav.targPickHistory==targAmpPick)))>=e.numAdjustTrials %if there are enough 'completed' trials, and the last trial was 'scorable' (this latter condition keeps it from running down repeatedly)
        trialSubset = behav.score(behav.targPickHistory==targAmpPick&~isnan(behav.score));
        score = mean(trialSubset(end-e.numAdjustTrials+1:end));
        %             fprintf('\nScore for adjustment: %f',score);
        if score<min(e.behavLimits) %if the score is below the accepted range
            if targAmpPick==1&&behav.targAmp(1)>=e.targAmpStep
                behav.targAmp(1)=behav.targAmp(1)-e.targAmpStep; %make the 'slow' target slower (down to zero)
            elseif targAmpPick==2
                behav.targAmp(2)=behav.targAmp(2)+e.targAmpStep; %make the fast target faster
            end;
        elseif score>max(e.behavLimits) %if the score is above the accepted range
            if targAmpPick==1&&behav.targAmp(1)<1-(2*e.targAmpStep)
                behav.targAmp(1)=behav.targAmp(1)+e.targAmpStep; %make the slow target faster
            elseif targAmpPick==2&&behav.targAmp(2)>=1+(2*e.targAmpStep)
                behav.targAmp(2)=behav.targAmp(2)-e.targAmpStep; %make the fast target slower (down to 1+e.targAmpStep)
            end;
        end;
    end;
end;
trialData{6} = sprintf('Hits: %d, Incorrects: %d (%.1f%%), Other: %d',sum(behav.score==1&~isnan(behav.score)),sum(behav.score==0&~isnan(behav.score)),nanmean(behav.score)*100,sum(isnan(behav.score)));
trialData{7} = sprintf('Target amplitudes: %.2f, %.2f',behav.targAmp);
%     trialData{7} = ['RTs (ms): ' sprintf('%4.1f, ',fliplr(behav.RT(~isnan(behav.RT))*1000))];
