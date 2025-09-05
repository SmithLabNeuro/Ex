function result = ex_SimpleTaskDemo(e)
% ex file: ex_SimpleTaskDemo
%
%
% This is a simple task file for a Memory-Guided Saccade - it's overly
% simplified and heavily commented to make it an instructive demo, but does
% work as a basic task.
%
% XML REQUIREMENTS - these are the parameters that are required in an XML
% file that points to this ex-function
%
% distance: the distance of the target in pixels from the fixation point
% angle: angle of the target dot from the fixation point 0-360 size: the
% size of the target in pixels fixX: X-location of the fixation point fixY:
% Y-location of the fixation point fixRad: radius of the fixation point in
% pixels targetColor: a 3 element [R G B] vector for the target color
% targetDuration: duration that target is on screen (ms) noFixTimeout:
% timeout punishment for aborted trial (ms) noChoiceTimeout: timeout
% punishement for failing to leave fixation (ms) preTargetFixation time
% required for fixation before the target flashes saccadeInitiate: maximum
% time allowed to leave fixation window saccadeTime: maximum time allowed
% to reach target stayOnTarget: length of target fixation required
% postTargetFixation: length of fixation required after target is off
% (before go cue of fixation offset)
%

%% Some initial setup

global params codes behav;
% params contains all the global variables that are set by the ex-function
% (both overall globals, as well as items that may be set in the rig or
% subject XML files)
%
% codes is the struct that has the names of special codes used for sending
% digital information to the data collection computers.
%
% behav is a place to store data that will persist across trials. It is
% empty by default and you need to fill it in with needed data

% take radius and angle and figure out x/y for saccade direction
theta = deg2rad(e.angle);
newX = round(e.distance*cos(theta));
newY = round(e.distance*sin(theta));

%% This part of the function has some communication with showex.

%    In this short block of code we setup the visual stimuli we want to
%    use. This is accomplished by 'msg' commands in which runex
%    communicates with showex. Here we're using 'set' commands, which tells
%    showex to setup (but not yet display) certain objects with parameters
%    indicated by the passed variables. These set commands make 'oval'
%    objects, but there are a variety of object types that showex can
%    handle (indicated by the stim_XXX.m functions - there is a
%    "stim_oval.m" function for example).

% This is the object number we want the diode "attached" to. Attaching the
% diode to an object means the diode will be flash white every time the
% object is turned on or off.
diodeObjID = 2; 

% Object 1 is the fixation point. It's important to note that objects are
% drawn in reverse order. So if you want the fixation point to always be
% "on top" of everything, make it object 1. That way if you happen to put
% another object in the same location, the fixation point will be "on top"
% of it.
msg('set 1 oval 0 %i %i %i %i %i %i',[e.fixX e.fixY e.fixRad e.fixColor(1) e.fixColor(2) e.fixColor(3)]);

% Object 2 is the target
msg('set 2 oval 0 %i %i %i %i %i %i',[newX newY e.size e.targetColor(1) e.targetColor(2) e.targetColor(3)]);

% This command tells showex to attached the diode to object 2, so it's
% flashed when object 2 is turned on or off
msg(['diode ' num2str(diodeObjID)]);

% Note that all the commands above are 'msg' commands, and we will use a
% 'msgAndWait' below. 'msg' commands are non-blocking - runex sends them
% and doesn't wait for showex to perform them.


%% Actual Trial Activities Begin Here

% This is a command that tells shows to turn on object 1 (the fixation
% point). Because this is a "msgAndWait" command, the function will not
% return until showex actually swaps the graphics buffer to turn on that
% stimulus. This, the return of this function is very precisely timed to
% the actual appearance of object 1 on the screen.
msgAndWait('obj_on 1');

% This code sends a digital code to the data collection system, and also
% stores that code in the behavioral ".mat" file on the runex computer.
% It's very important that this code be sent *after* the msgAndWait command
% in this instances, because it gives the best alignment of that code in
% the data with the actual appearance of the visual stimulus.
sendCode(codes.FIX_ON);

% Now that the fixation spot is on the screen, we need to wait for the
% subject to fixate. We use "waitForFixation" to check for that. However,
% this if statement is formulated with a "not" in front. So we should read
% this statement as "if the subject doesn't fixate, we end the trial".
if ~waitForFixation(e.timeToFix,e.fixX,e.fixY,params.fixWinRad);
    % If the subject failed to achieve fixation
    sendCode(codes.IGNORED);
    msgAndWait('all_off'); % this turns off all graphics on the screen (all object numbers)
    sendCode(codes.FIX_OFF);
    waitForMS(e.noFixTimeout); % this is essentially the same as a "pause"
    result = codes.IGNORED;
    return;
end
sendCode(codes.FIXATE); % at this point the eyes have entered the fixation window

% Now that the subject has reached the fixation window, we should wait for
% a period of time until we turn on the target flash.
if ~waitForMS(e.preTargetFixation,e.fixX,e.fixY,params.fixWinRad)
    % hold fixation before stimulus comes on
    sendCode(codes.BROKE_FIX);
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    waitForMS(e.noFixTimeout);
    result = codes.BROKE_FIX;
    return;
end

% turn on the target and send a code after it's on
msgAndWait('obj_on 2');
sendCode(codes.TARG_ON);

% this is to keep the target on for a specified duration of time. And
% important note here is that because of latencies between runex and
% showex, you might ask for a specific amount of time but the target will
% certainly be on for a slightly longer period of time (usually a frame or
% two). If you want the target on for a very specific amount of time, you
% need to set the "frames" parameter when you create the object, and you
% can use the "waitForDisplay" function to let showex count the time the
% target stays on. That provides more precise object timing.
if ~waitForMS(e.targetDuration,e.fixX,e.fixY,params.fixWinRad)
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

% This is an additional fixation period required after the target ends. In
% this task, this is the "memory delay", because the subject will need to
% make a saccade to the remembered target location
if ~waitForMS(e.postTargetFixation,e.fixX,e.fixY,params.fixWinRad)
    sendCode(codes.BROKE_FIX);
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    waitForMS(e.noFixTimeout);
    result = codes.BROKE_FIX;
    return;
end

% Now, turn off the fixation point. This is the cue to let the subject make
% a saccade to the remembered target location
msgAndWait('obj_off 1');
sendCode(codes.FIX_OFF);

% This is a special parameter we can enable. If 'recenterFixWin' is
% enabled, then the current eye position at the function on set will be
% used to monitor the position of the eye to determine a saccade, instead
% of using the actual fixation dot location. This is important to help be
% sure to detect small saccades and handle any fixational drift.
if params.recenterFixWin
    newFixWinRad = params.sacWinRad;
else
    newFixWinRad = params.fixWinRad;
end

saccadeInitiateRemainder = e.saccadeInitiate;

% This is a nice way to detect an optional parameter in a XML file. If it's
% not a field of the struct, the parameter didn't exist and therefore the
% code should handle that gracefully. This is often a good choice of how to
% handle backwards compatibility if you want to add a new parameter to a
% XML file but still have the code work with older XML files that didn't
% specify that parameter.
if isfield(e,'minRT') 
    
    % if requested, wait for a minimum amount before saccade (because very
    % short reaction times meant the subject didn't perceive the go cue but
    % instead "guessed" when it would go away). Here, we use the
    % "waitForMS" function because we want the subject to keep their eyes
    % on the fixation dot during this time, and the function returns if
    % they *leave* the fixation window. So this if/then statement can be
    % read as "if the subject didn't stay fixating for minRT, end the
    % trial".
    if ~waitForMS(e.minRT,e.fixX,e.fixY,params.fixWinRad)
        sendCode(codes.SACCADE);
        sendCode(codes.BROKE_FIX);
        msgAndWait('all_off');
        waitForMS(e.noFixTimeout);
        result = codes.BROKE_FIX;
        return;
    end

    % calculate how much time left to wait
    saccadeInitiateRemainder = saccadeInitiateRemainder-e.minRT;
end

% Again, use a waitForMS, but in a different way. Here, we wait for the
% remainder of the saccadeInitiate time, and if the subject stayed in the
% window the whole time, that means they did not make a saccade. That's an
% error in this task ("NO_CHOICE") and means the trial should end.
if waitForMS(saccadeInitiateRemainder,e.fixX,e.fixY,newFixWinRad,'recenterFlag',params.recenterFixWin)
    sendCode(codes.NO_CHOICE);
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    waitForMS(e.noChoiceTimeout); % Here we have a specific pause related to this condition, which we call a timeout
    result = codes.NO_CHOICE;
    return;
end

% If the subject got to here, it means they left the fixation window before
% the timer elapsed. So, we indicate this was a saccade. Of course offline
% analysis is necessary to determine if it's truly a saccade vs. some other
% accidental move of the eyes.
sendCode(codes.SACCADE);

% Now, if you left the fixation window we need the subject to get to the
% target window within a specific amount of time.
if ~waitForFixation(e.saccadeTime,newX,newY,params.targWinRad)
    % didn't reach target
    sendCode(codes.NO_CHOICE);
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    waitForMS(e.noChoiceTimeout); % timeout
    result = codes.NO_CHOICE;
    return;
end

% this code means the subject reached the target window
sendCode(codes.ACQUIRE_TARG);

% here we require the subject to stay in the target window for a bit of
% time so that we don't let them just barely skim through it for a brief
% period of time. They have to stop and hold for a short period.
if ~waitForMS(e.stayOnTarget,newX,newY,params.targWinRad)
    sendCode(codes.BROKE_TARG);
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    result = codes.BROKE_TARG;
    return;
end

% If we get here, the trial's a success. So, let's send a bunch of
% additional codes
sendCode(codes.FIXATE);
sendCode(codes.CORRECT);
sendCode(codes.REWARD);

% Go ahead and reward the subject at this point - success!
giveJuice();
result = codes.CORRECT;

% It may be convenient to have a little bit of time between trials just to
% keep things from going too fast. So this is an optional parameter for
% that. Since it's after the reward, the subject doesn't know that this is
% happening at the end of this trial vs. at the beginning of the next -
% there's nothing on the screen.
if isfield(e,'InterTrialPause')
    waitForMS(e.InterTrialPause);
end

