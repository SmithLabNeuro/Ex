function result = ex_SimpleJoystickTaskDemo(e)
% ex file: ex_SimpleJoystickTaskDemo
%
%
% This is a simple task file for a Center-Out Joystick Reach - it's overly
% simplified and heavily commented to make it an instructive demo, but does
% work as a basic task.
%
% XML REQUIREMENTS - these are the parameters that are required in an XML
% file that points to this ex-function
%
% distance: the distance of the target in pixels from the fixation point
% angle: angle of the target dot from the fixation point 0-360 size: the
% size of the target in pixels 
% fixX: X-location of the fixation point 
% fixY: Y-location of the fixation point 
% fixRad: radius of the fixation point in pixels 
% cursorRad: radius of the cursor in pixels
% targetColor: a 3 element [R G B] vector for the target color
% cursorColor: a 3 element [R G B] vector for the cursor color
% targetDuration: duration that target is on screen (ms) 
% noFixTimeout: timeout punishment for aborted trial (ms) 
% noChoiceTimeout: timeout punishement for failing to leave fixation (ms) 
% preTargetFixation time required for fixation before the target flashes 
% joystickInitiate: maximum time allowed to leave fixation window 
% joystickTime: maximum time allowed to reach target 
% stayOnTarget: length of target fixation required (ms)
% postTargetFixation: length of fixation required after target is off
% (ms, before go cue of fixation offset)
%
% Default joystick parameters in the XML:
% joystickWinRad: radius of the joystick acquisition window in pixels
% screenPixelLimit: maximum pixel value for the joystick

%% Some initial setup

global params codes;
% params contains all the global variables that are set by the ex-function
% (both overall globals, as well as items that may be set in the rig or
% subject XML files)
%
% codes is the struct that has the names of special codes used for sending
% digital information to the data collection computers.

%% This part of the function sets up event checking and object parameters

% first for the reach target: 
% take radius and angle and figure out x/y for saccade direction
theta = deg2rad(e.angle);
newX = round(e.distance*cos(theta));
newY = round(e.distance*sin(theta));

% next for the cursor:
cursorObjID = 3; % this is the object number for the cursor
cursorColor = e.cursorColor; % this is the color of the cursor

% set up event checker requirements (collect function parameters into cell arrays for use in waitForEvent)

% 1) fixation requirements
fixationArgs = {e.fixX, e.fixY, params.fixWinRad};

% 2) Hall Effect Joystick requirements
% for acquiring fixation, we want the joystick to be in the fixation window
joystickHEAcqPosArgs = {e.fixX, e.fixY, e.joystickWinRad, e.screenPixelLimit}; 
% for holding fixation, we want the joystick to be in the fixation window
joystickHEHoldMsArgs = [joystickHEAcqPosArgs, cursorObjID, e.cursorRad, cursorColor, false]; % the last argument is a flag to draw the cursor or not
% for moving the joystick from the initial fixation position, also uses the cursor parameters to draw the cursor on the screen
joystickHEMoveArgs = {[e.fixX, e.fixY], e.joystickWinRad, e.screenPixelLimit, cursorObjId, e.cursorRad, cursorColor};
% for moving the joystick to the target, also uses the cursor parameters to draw the cursor on the screen
joystickHETargetArgs = {newX, newY, cursorObjID, e.cursorRad, cursorColor, e.targWinRad, e.screenPixelLimit};
% for holding the joystick in the target window, also uses the cursor parameters to draw the cursor on the screen
joystickHEHoldTargetArgs = {newX, newY, e.targWinRad, e.screenPixelLimit, cursorObjID, e.cursorRad, cursorColor, true, e.stayOnTarget};

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
diodeObjID = cursorObjID; 

% Object 1 is the fixation point. It's important to note that objects are
% drawn in reverse order. So if you want the fixation point to always be
% "on top" of everything, make it object 1. That way if you happen to put
% another object in the same location, the fixation point will be "on top"
% of it.
msg('set 1 oval 0 %i %i %i %i %i %i',[e.fixX e.fixY e.fixRad e.fixColor(1) e.fixColor(2) e.fixColor(3)]);

% Object 2 is the target
msg('set 2 oval 0 %i %i %i %i %i %i',[newX newY e.size e.targetColor(1) e.targetColor(2) e.targetColor(3)]);

% Object 3 is the cursor
msg('set 3 oval 0 %i %i %i %i %i %i',[e.fixX e.fixY e.cursorRad e.cursorColor(1) e.cursorColor(2) e.cursorColor(3)]);

% This command tells showex to attached the diode to object 2, so it's
% flashed when object 2 is turned on or off
msg(['diode ' num2str(diodeObjID)]);
msgAndWait('ack');

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
% subject to fixate. Fixation in the joystick task entails maintiaining eye 
% position at the fixation point, and also moving the joystick to the fixation point. 
% To check multiple conditions, we use waitForEvent.
% waitForEvent takes a time limit, a cell array of function handles, and a cell array of arguments for those functions. 
% It will return true when all of the functions return true within the time limit. 
% In this case, we are checking for either joystick acquisition and fixation acquisition.
[acquiredFixation, ~] = waitForEvent(e.timeToFix, {@joystickAcquire, @fixationAcquire}, {joystickHEAcqPosArgs, fixationArgs});
if ~acquiredFixation
    % If the subject failed to achieve fixation
    result = codes.IGNORED;
    sendCode(result);
    msgAndWait('all_off'); % this turns off all graphics on the screen (all object numbers)
    sendCode(codes.FIX_OFF);
    waitForMS(e.noFixTimeout); % this is essentially the same as a "pause"
    return
end
sendCode(codes.FIXATE); % at this point the eyes have entered the fixation window

% Now that the subject has reached the fixation window, we should wait for
% a period of time until we turn on the target flash.
% Now we use waitForEvent again, but this time we want to hold fixation for a period of time, so we pass in new event checker function handles.
joystickHEHoldMsArgs{end+1} = e.preTargetFixation; % add the time to hold fixation before target flash
[heldFixation, ~] = waitForEvent(e.preTargetFixation, {@joystickHold, @fixationHold}, {joystickHEHoldMsArgs, fixationArgs});
if ~heldFixation
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

% this is to keep the target on for a specified duration of time. 
% since the behavior doesn't change (still fixating), the same event checker functions are used 
joystickHEHoldMsArgs{end} = e.targetDuration; % update the time to hold fixation during target flash
[heldFixation, ~] = waitForEvent(e.targetDuration, {@joystickHold, @fixationHold}, {joystickHEHoldMsArgs, fixationArgs});
if ~heldFixation
    sendCode(codes.BROKE_FIX);
    msgAndWait('all_off');
    sendCode(codes.TARG_OFF);
    sendCode(codes.FIX_OFF);
    waitForMS(e.noFixTimeout);
    result = codes.BROKE_FIX;
    return;
end

% Now, turn off the fixation point. This is the cue to let the subject make
% a saccade to the remembered target location
msgAndWait('obj_off 1');
sendCode(codes.FIX_OFF);

% Again, use a waitForMS, but in a different way. Here, we wait for the
% remainder of the saccadeInitiate time, and if the subject stayed in the
% window the whole time, that means they did not make a saccade. That's an
% error in this task ("NO_CHOICE") and means the trial should end.
[reacted, ~] = waitForEvent(e.joystickInitiate, {@joystickMove}, {joystickHEMoveArgs});
if ~reacted
    sendCode(codes.NO_CHOICE);
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    waitForMS(e.noChoiceTimeout); % Here we have a specific pause related to this condition, which we call a timeout
    result = codes.NO_CHOICE;
    return;
end

% If the subject got to here, it means they left the fixation window before
% the timer elapsed. Even though it is technically a joystick movement, we will call 
% it a saccade for the purposes of this task. So, we send a code to indicate that the subject made a saccade.
sendCode(codes.SACCADE);

% Now, if you left the fixation window we need the subject to get to the
% target window within a specific amount of time.
[reachTarget, ~] = waitForEvent(e.joystickTime, {@joystickMoveToTarget}, {joystickHETargetArgs});
if ~reachTarget
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
[heldTarget, ~] = waitForEvent(e.stayOnTarget, {@joystickHold}, {joystickHEHoldTargetArgs});
if ~heldTarget
    sendCode(codes.BROKE_TARG);
    msgAndWait('all_off');
    sendCode(codes.FIX_OFF);
    result = codes.BROKE_TARG;
    return;
end

% If we get here, the trial's a success. So, let's send a bunch of
% additional codes
sendCode(codes.FIXATE); % this is a bit of a misnomer, but it means that fixation was met at the target
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

