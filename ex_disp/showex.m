function showex(varargin)
% function showex('param1','value1',...)
%
% Main function for displaying stimuli for Ex.
%
% Parameters:
% 'ip' (string): the IP address for the control machine.
% 'photometer' (string): the file name for the photometer values text file
% used for gamma correction.
%
% Modified 16Sep2013 - (1) auto detect machine name and use that to pick the
% individualized photometer files. (2) change inputs to parameter-value
% pairs and set defautls. Added the photometer text file input argument
% option. (3) added 'eval_str' command, which will evaluate strings passed
% from the control side. This is used to add local directories, etc. (4)
% added directory handling stuff to make sure the required directories are
% on the search path. (5) put in some cells, comments, etc. -ACS
%
% Modified 04Sep2013 - 'dset' flashes the diode on a set command. -ACS
%
% Modified 17Jul2013 - made diode slightly smaller to prevent light leaking
% out the side. -ACS
%
% Modified 14Dec2015 - Using KbQueue commands instead of GetChar/ListenChar
%
% Modified 24Feb2016 - added obj_switch and the ability to toggle more than
% one object on or off in a single command. Also added Priority features
% for realtime operation, and Verbosity option to suppress display. 
% General cleanup as well. -MAS
%

%% Start a keyboard queue
KbQueueCreate;
KbQueueStart;

%% initialize global variables
global objects;
global sv;
global w;

%% Startup message
thisFileInfo = dir(which(mfilename));
fprintf('*****Running SHOWEX (last modified %s). Started %s*****\n\n',thisFileInfo.date,datestr(now));

%% Figure out our machine name
[~,sv.machine] = system('hostname'); %name of the showex machine -ACS 16Sep2013
sv.machine = lower(deblank(cell2mat(regexp(sv.machine,'^[^\.]+','match'))));
assert(~isempty(sv.machine),'Machine name must not be empty'); 
assert(sum(isspace(sv.machine))==0,'Machine name must not have spaces');

%% input stuff
p = inputParser;
p.addParameter('ip','192.168.1.11',@ischar);
p.addParameter('photometer',sprintf('photometer_%s.txt',sv.machine),@ischar);

p.parse(varargin{:});

ipAddress = p.Results.ip;
photometerFile = p.Results.photometer;

%% directory stuff:
thisFile = mfilename('fullpath');
[runexDirectory,~,~] = fileparts(thisFile); clear thisFile;
runexDirectory = regexp(runexDirectory,sprintf('.*(?!\\%c).*\\%c',filesep,filesep),'match'); %everything up to and including the last file separator
requiredDirectories = {'ex_disp','stim','ex_control'};
for dx = 1:numel(requiredDirectories)
    addpath([runexDirectory{:},requiredDirectories{dx}]);
end

%% local directory
exGlobals;
if exist([params.localExDir,filesep,'display'],'dir')
    addpath([params.localExDir,filesep,'display']);
end

%% UDP initialization
delete(instrfind)
%u = udp(ipAddress,'RemotePort',8844,'LocalPort',8866);
%fopen(u);
%local_ip = char(java.net.InetAddress.getLocalHost.getHostAddress);
matlabUDP2('all_close');
socket = matlabUDP2('open','192.168.1.10', ipAddress, 4243);
matlabUDP2('send',socket,'ack');

%% Put our code in high-priority/realtime mode
origPriority = Priority(1); % 1 or > is realtime mode; 0 is normal
newPriority = Priority(1); % check that the priority is now 1
if newPriority ~= 1 
        ansr = questdlg('Unable to change Priority. Probably PsychLinuxConfiguration was not run. Continue?');
    switch lower(ansr(1))
        case 'y'
            warning('Running without setting Priority');
        otherwise
            return; %punt
    end
end

%% Set verbosity. 3 is default, 1 means only output errors
origVerbosity = Screen('Preference', 'Verbosity', 1);

%% screen initialization
Screen('Preference', 'SkipSyncTests', 0);

% check to see if it can load photometerFile
% if not, gui dialog to ask if you want to load the default file
if exist(photometerFile,'file')
    gam = load(photometerFile);
else
    ansr = questdlg('No local photometer file was found, do you want to load defaults? (if not, press no or cancel and find the file)');
    switch lower(ansr(1))
        case 'y'
            photometerFile = 'photometervals.txt';        
            gam = load(photometerFile);
        otherwise
            return; %punt
    end
end

vals = makeGammaTable(gam(:,1),gam(:,2));
Screen('LoadNormalizedGammaTable', 0, vals);

%% Other initialization
% specifies size and location of photodiode square:
% upper left of screen, 30 x 30 pixels
diodeLoc = [10 10 40 40];

AssertOpenGL;

screens = Screen('Screens');
screenNumber = max(screens);

sv.white = WhiteIndex(screenNumber);
sv.black = BlackIndex(screenNumber);
sv.gray = (sv.white+sv.black)/2;
if round(sv.gray) == sv.white
    sv.gray = sv.black;
end
sv.inc = sv.white-sv.gray;
sv.bgColor = sv.gray;

% open the display window
[w, sv.screenRect] = Screen('OpenWindow',screenNumber,sv.bgColor);
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

sv.midScreen = sv.screenRect(3:4)/2;

% show the photodiode square
Screen(w,'FillRect',sv.black,diodeLoc);
Screen('Flip',w);

%% get the frame time (inter-frame interval)
sv.ifi = Screen('GetFlipInterval',w);
tvblt = inf;
timing_flag = false;
diode_timing_enable = false;
op = [];
%save('/home/smithlab/Dropbox/smithlabrig/testData/21.mat','op');

%% turn off mouse pointer
HideCursor();

%% pre-allocate object-related variables
maxObjects = 20;
objects = cell(maxObjects,1);
visible = zeros(maxObjects,1);
visibleQueue = zeros(maxObjects,1);

%% intialize some flags/values for the display loop
queueing = false;
diodeVal = 0;
diodeFlashFull = 255;
diodeFlashPartial = 192;
diodeValToFlash = diodeFlashPartial;
diodeObj = 1;
returnReceiptPending = false;
doneMsgPending = false;
doneMsg = '';

%% DISPLAY LOOP
while true % display loop
    try % master try statement to catch errors in the display loop
        
        if matlabUDP2('check',socket) % received network message
            [s1, s] = strtok(matlabUDP2('receive',socket));
            returnReceiptPending = true;
            
            %for debugging use this line
            %disp([s1 s]);
            
            switch s1
                case {'set','dset'}
                    % arguments: object id + parameters for the stim_ function
                    [objID, withoutID] = strtok(s);
                    objID = str2double(objID);
                    [objType, args] = strtok(withoutID);
                    if strcmp(s1,'dset')
                        diodeVal=0;
                        feval(['stim_',objType],'setup',w,objID,args);                        
                    else
                        feval(['stim_',objType],'setup',w,objID,args);
                    end
                    
                case 'mset'
                    % arguments: object id + parameters for the stim_ function                  
                    diodeVal=0;
                    sSplit = strsplit(s,'set');
                    nset = length(sSplit);
                    
                    for I=1:nset
                       [objID, withoutID] = strtok(sSplit{I});
                       [objType, args] = strtok(withoutID);
                       objID = str2double(objID);
                       feval(['stim_',objType],'setup',w,objID,args);
                    end
                    
                case 'obj_on'
                    % arguments: (1-n) objects ids to make visible
                    args = cell2mat(textscan(s,'%n'));
                    if queueing
                        visibleQueue(args) = 1;
                    else
                        visible(args) = 1;
                    end
                    
                case 'obj_off'
                    % arguments: (1-n) objects ids to turn off
                    args = cell2mat(textscan(s,'%n'));
                    if queueing
                        visibleQueue(args) = 0;
                    else
                        visible(args) = 0;
                    end
                    
                case 'obj_switch'
                    % arguments: (1-n) objects ids to turn on (if positive) or off (if negative)
                    args = cell2mat(textscan(s,'%n'));
                    if queueing
                        visibleQueue(args(args>0)) = 1;
                        visibleQueue(-args(args<0)) = 0;
                    else
                        visible(args(args>0)) = 1;
                        visible(-args(args<0)) = 0;
                    end
                    
                case 'all_on'
                    % arguments: none
                    if queueing
                        visibleQueue(~cellfun(@isempty,objects)) = 1;
                    else
                        visible(~cellfun(@isempty,objects)) = 1;
                    end
                    
                case 'rem_all'
                    % arguments: none
                    for objID = 1:length(objects)
                        if ~isempty(objects{objID})&&isfield(objects{objID},'type')
                            feval(['stim_',objects{objID}.type],'cleanup',w,objID);
                        end
                    end
                    objects = cell(maxObjects,1);
                    visible(:) = 0;
                    queueing = false;

                case 'all_off'
                    % arguments: none
                    if queueing
                        visibleQueue(:) = 0;
                    else
                        visible(:) = 0;
                    end
                case 'timing_begin'
                    % arguments: none
                    %timing_flag = true;
                    dropCount = 0;
                    frameCount = 0;
                    frameDrop = 0;
                    
                case 'timing_end'
                    % arguments: none
                    %timing_flag = false;
                    % send number of drops here, and reset it to zero after
                    s1 = num2str(dropCount);

                case 'queue_begin'
                    % arguments: none
                    queueing = true;
                    visibleQueue = visible;
                    
                case 'queue_end'
                    % arguments: none
                    queueing = false;
                    visible = visibleQueue;
                    
                case 'bg_color'
                    % arguments: 3 numbers indicating the background color
                    args = textscan(s,'%n');
                    sv.bgColor = args{1};
                    Screen(w,'FillRect',sv.bgColor);
                    
                case 'diode'
                    % arguments: object ID to tie the diode square to
                    % diodeObj = str2double(s);
                    diodeObj = cell2mat(textscan(s,'%n'));
                    
                case 'diode_setup'
                    % arguments: left/top/right/bottom of diode position
                    args = textscan(s,'%n');
                    diodeLoc = [args{1}(1) args{1}(2) args{1}(3) args{1}(4)];
                    Screen(w,'FillRect',sv.black,diodeLoc);
                    
                case 'diode_timing'
                    % set the flag here
                    diode_timing_enable = true;
                    % maybe just an option on diode
                    
                case 'framerate'
                    % arguments: none
                    % because we re-set s1 here, this value gets
                    % sent back instead of the word 'framerate'
                    s1 = sprintf('%.6f', Screen('GetFlipInterval',w));
                    
                case 'resolution'
                    % arguments: none
                    % because we re-set s1 here, this value gets
                    % sent back instead of the word 'resolution'
                    res = Screen('Resolution',w);
                    s1 = [num2str(res.width),' ',num2str(res.height),' ',num2str(res.pixelSize),' ',num2str(res.hz)];
                    
                case 'screen'
                    % arguments: screen distance and pixpercm
                    args = textscan(s,'%n');
                    sv.scrd = args{1}(1);       % screen distance in cm
                    sv.pixpercm = args{1}(2);   % pixels per cm
                    sv.ppd = tan(deg2rad(1)) * sv.scrd * sv.pixpercm; % pixels per degree
                    
                case 'eval_str' %added to evaluate functions in showex remotely (e.g. addpath)
                    eval(s);
                    
            end
            
            % Now we've gotten our command, move on to displaying
            % the objects but don't forget to send back the return receipt
        end
        
        vis = find(visible);
        
        % display section
        if ~isempty(vis)
            for i = length(vis):-1:1
                objID = vis(i);
                try % call the stim function
                    feval(['stim_',objects{objID}.type],'display',w,objID); 
                    % only timed objects have fc set to a non-zero value
                    if objects{objID}.frame+1 == objects{objID}.fc
                        visible(objID) = 0;
                        if queueing, visibleQueue(objID) = 0; end
                        doneMsgPending = true;
                        doneMsg = sprintf('done %i', objID);
                        objects{objID}.frame = 0;
                        diodeValToFlash = diodeFlashFull;
                    else
                        objects{objID}.frame = objects{objID}.frame+1;
                    end
                catch ME %for debugging -ACS 11Jun2015
                    switch ME.identifier
                      case 'MATLAB:nonStrucReference' % objects{objID} is empty -ACS 11Jun2015
                        warning('SHOWEX:callToDeletedObject','An attempt was made to draw or clear a deleted object... resuming...');
                        continue;
                      otherwise
                        rethrow(ME);
                    end
                end
            end
        end %closes if ~isempty(vis) condition
        
        % Make diode switch on alternate frames:
        % White on 1st frame when object is on, then alternate between black and 
        % mid-gray during the stimulus, then white on 1st frame when object is off
        % if visible(diodeObj)
        if any(visible(diodeObj))
            switch diodeVal
              case 0
                diodeVal = diodeFlashFull;
                if diode_timing_enable
                    timing_flag = true;
                end
              case diodeFlashFull
                diodeVal = 1;
              case diodeFlashPartial
                diodeVal = 1;
              otherwise
                diodeVal = mod(diodeVal + 128,256);
            end
        else
            switch diodeVal
              case diodeFlashFull
                diodeVal = 0;
              case diodeFlashPartial
                diodeVal = 0;
              case 0
                diodeVal = 0;
              otherwise
                diodeVal = diodeValToFlash;
            end
            timing_flag = false;
        end
        
        Screen(w,'FillRect',diodeVal * [1 1 1],diodeLoc);
        
%       Screen('Flip',w); % Main Flip Command
        [tvbl,~,~,~,~] = Screen('Flip',w); % Main Flip Command
       
        if timing_flag
            frameCount = frameCount+1;
            if tvbl - tvblt > sv.ifi * 1.5
                %disp(tvbl-tvblt);
                %timing_flag = false;
                tvblt = Inf; 
                dropCount = dropCount + 1;
                frameDrop(dropCount) = frameCount;
                fprintf('dropping %d th frame',frameCount);
            end
        end
       tvblt = tvbl;
        % Right after the Flip, send back the return message FIRST
        if returnReceiptPending 
            matlabUDP2('send', socket,s1);
            returnReceiptPending = false;
        end
        
        % Right after the Flip, send back the done message
        if doneMsgPending
            matlabUDP2('send',socket, doneMsg);
            doneMsg = '';
            doneMsgPending = false;
            diodeValToFlash = diodeFlashPartial;
        end
        
        % check for keyboard input
        [ keyIsDown, keyCode] = KbQueueCheck;
        if keyIsDown
            c = KbName(keyCode);
            KbQueueFlush;
            if numel(c)>1; continue; end %keyboard mash and other weirdness
            
            if c == 'x'
                break;
            end
        end
        
    catch ERR
        matlabUDP2('send', socket,'abort'); % problem with UDP communication
        disp(ERR.message);
        for sx = 1:numel(ERR.stack)
            disp(ERR.stack(sx));
        end
        Screen(w,'FillRect',sv.bgColor);
        Screen('Flip',w);
        fprintf('Waiting to resume... %s\n',datestr(now));
        tic
        while true
            [ keyIsDown, keyCode] = KbQueueCheck;
            if keyIsDown
                c = KbName(keyCode);
                KbQueueFlush;
                if numel(c)>1; continue; end %keyboard mash and other weirdness
                
                if c == 'x'
                    break;
                end
            end
            
            if toc > 60
                sca
                error('Timeout - waited for 1 minute with no resume command');
            end
            
            if matlabUDP2('check',socket) % received network message
                [s1, s] = strtok(matlabUDP2('receive',socket));
                
                %for debugging use this line
                disp([s1 s]);
                switch s1
                    case 'resume'
                        %clear everything and return to the display loop:
                        objects = cell(maxObjects,1);
                        visible(:) = 0;
                        queueing = false;
                        returnReceiptPending = false;
                        doneMsgPending = false;
                        doneMsg = '';
                        fprintf('Resuming... %s\n',datestr(now));
                        matlabUDP2('send', socket,'resume'); %handshake with runex
                        break;
                    otherwise
                        %keep waiting
                end
            end
        end
    end
end

%% Clean up on exit:
Priority(origPriority);
Screen('Preference', 'Verbosity', origVerbosity);
sca
matlabUDP2('all_close');

end

% function cleanUpObj(o)
% switch o.type
%     case 10
%         Screen('Close',o.grating);
%         Screen('Close',o.mask);
% end
% end
