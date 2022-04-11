function [success, funcSuccess, extraFuncOutput] = waitForEvent(timeout, eventCheckerFunction, eventFunctionInputs, successIfTimeElapsed)
% function waitForEvent(timeout)
%
% function to generalize waitFor functionality:
%
% timeout: how long to wait for the event to happen
%
% eventCheckerFunction: cell of function handles with API
%   [success] = eventCheckerFunction{i}(eventFunctionInputs{i}{:})
% so that it returns a success signal through every loop. If you want the
% function to retain some history through the loop, one solution is to use
% persistent variables in the eventCheckerFunction
% NOTE: eventCheckerFunction can also be a single anonymous function, at
% which point eventFunctionInputs should be just a cell array of the inputs
% to that function (see below), and the line above would not have the index
% i:
%   [success] = eventCheckerFunction(loopStart, loopTop, eventFunctionInputs{:})
%
% eventFunctionInputs: cell of cells of function inputs to each anonymous
% function eventCheckerFunction. So for two functions with three arguments
% each, it would look like 
%   {{fun1arg1, fun1arg2, fun1arg3},{fun2arg1,fun2arg2,fun2arg3}}
% NOTE: if eventCheckerFunction is instead a single anonymous function (not
% a cell), then eventFunctionInputs should just be a single cell with the
% inputs of eventCheckerFunction:
%   e.g. {fun1arg1, fun1arg2, fun1arg3}

global debug params sockets;

if nargin<1, timeout = 10; end;


if ~iscell(eventCheckerFunction)
    eventCheckerFunction = {eventCheckerFunction};
    eventFunctionInputs = {eventFunctionInputs};
end

loopStart = GetSecs;
ticStart = tic;
success = 0;
runOnceMore = true;
firstLoop = true;
while (toc(ticStart)*1000) <= timeout
    loopTop = GetSecs;
    
    funcSuccess = zeros(1,length(eventCheckerFunction));
    displayUpdateString = cell(1,length(eventCheckerFunction));
    fixWindowUpdateInputs = cell(1, length(eventCheckerFunction));
    extraFuncOutput = cell(1, length(eventCheckerFunction));
    for funcInd = 1:length(eventCheckerFunction)
        if nargout(eventCheckerFunction{funcInd}) == 4
            [funcSuccess(funcInd), displayUpdateString{funcInd}, fixWindowUpdateInputs{funcInd}, extraFuncOutput{funcInd}] = eventCheckerFunction{funcInd}(loopStart, loopTop, eventFunctionInputs{funcInd}{:});
        elseif nargout(eventCheckerFunction{funcInd}) == 3
            [funcSuccess(funcInd), displayUpdateString{funcInd}, fixWindowUpdateInputs{funcInd}] = eventCheckerFunction{funcInd}(loopStart, loopTop, eventFunctionInputs{funcInd}{:});
        elseif nargout(eventCheckerFunction{funcInd}) == 2
            [funcSuccess(funcInd), displayUpdateString{funcInd}] = eventCheckerFunction{funcInd}(loopStart, loopTop, eventFunctionInputs{funcInd}{:});
%             warning('%s needs a fixWindowUpdateInputs output (can be empty vector)', char(eventCheckerFunction{funcInd}));
        else
            funcSuccess(funcInd) = eventCheckerFunction{funcInd}(loopStart, loopTop, eventFunctionInputs{funcInd}{:});
            displayUpdateString{funcInd} = '';
%             warning('%s needs a displayUpdateString output (can be empty string)', char(eventCheckerFunction{funcInd}));
        end
    end
    
   
    % let's one press e.g. 'q' in order to quit a stimulus in the middle
    % with runex
    if keyboardEvents()
        success = 0;
        runOnceMore = false;
        break;
    end
    
    % ALL COMMANDS IN disaplyUpdateString MUST START WITH set
    for dispStrInd = length(displayUpdateString):-1:1
        if isempty(displayUpdateString{dispStrInd})
            displayUpdateString(dispStrInd) = [];
        elseif ~strcmp(displayUpdateString{dispStrInd}(1:3), 'set')
            error('bad display string %s:\n\nall display strings must start with set', displayUpdateString{dispStrInd})
        else
            % maybe add a white space at the end? might not matter for
            % these strings...
        end
    end
    if ~isempty(displayUpdateString)
%         disp(qx'hr')
        fullDisplayStr = strcat(displayUpdateString{:});
%         try
% fullDisplayStr
%     thisStart = tic;
if ~firstLoop
%     timeoutH=1;
%     warnOn = false;
%     while true
%         loopTopH = GetSecs;
% 
         while ~matlabUDP2('check',sockets(1))
         end
         if matlabUDP2('check', sockets(1))
            s = matlabUDP2('receive',sockets(1));
         end
%             toc(thisStart)
%             if debug
%                 fprintf('Rcvd: %s\n', s);
%             end;
%             if strcmp(s,'abort'), %added to potentially help with hangs. -ACS 03Sep2013
%                 error('waitFor:aborted','Abort signal received');
%             end;
%             if warnOn
%                 s={s, 'warn'};
%             end
%             break; %message received
%         elseif toc(thisStart)>timeoutH
%             error('waitFor:aborted','Timeout waiting for return message');
%         else
%             s = '';
%         end;
%         chk = GetSecs-loopTopH;
%         if (chk)>params.waitForTolerance, warning('waitFor:tooSlow','waitFor exceeded latency tolerance - %s - by %0.01f msecs',datestr(now), 1000*(chk - params.waitForTolerance)); warnOn = true;  end; %warn tolerance exceeded -acs22dec2012
%     end
else
    firstLoop = false;
end
          msg(['m' fullDisplayStr])
%           toc(thisStart)
%         catch err
%             a = 5;
%             b = 3;
%         end
    end
    
    fixXPerFunc = cell(1, length(fixWindowUpdateInputs));
    fixYPerFunc = cell(1, length(fixWindowUpdateInputs));
    sizeInfoPerFunc = cell(1, length(fixWindowUpdateInputs));
    maxSizeInput = 3;
    winColorsPerFunc = cell(1, length(fixWindowUpdateInputs));
    for fixWinUpdateInd = length(fixWindowUpdateInputs):-1:1
        if ~isempty(fixWindowUpdateInputs{fixWinUpdateInd})
            fixXPerFunc{fixWinUpdateInd} = fixWindowUpdateInputs{fixWinUpdateInd}{1};
            fixYPerFunc{fixWinUpdateInd} = fixWindowUpdateInputs{fixWinUpdateInd}{2};
            sizeInfoEvtFunc = fixWindowUpdateInputs{fixWinUpdateInd}{3};
            sizeInfoPerFunc{fixWinUpdateInd} = padarray(fixWindowUpdateInputs{fixWinUpdateInd}{3}, [maxSizeInput - size(sizeInfoEvtFunc, 1), 0], NaN, 'post');
            numWindows = length(fixXPerFunc{fixWinUpdateInd});
            winColorsEvtFunc = fixWindowUpdateInputs{fixWinUpdateInd}{4};
            if size(winColorsEvtFunc, 1) ~= numWindows
                winColorsEvtFunc = repmat(winColorsEvtFunc, numWindows, 1);
            end
            winColorsPerFunc{fixWinUpdateInd} = winColorsEvtFunc;
        else
            fixXPerFunc(fixWinUpdateInd) = [];
            fixYPerFunc(fixWinUpdateInd) = [];
            sizeInfoPerFunc(fixWinUpdateInd) = [];
            winColorsPerFunc(fixWinUpdateInd) = [];
        end
    end
    
    fixX = cat(2, fixXPerFunc{:});
    fixY = cat(2, fixYPerFunc{:});
    sizeInfo = cat(2, sizeInfoPerFunc{:});
    winColors = cat(1, winColorsPerFunc{:});
    drawFixationWindows(fixX, fixY, sizeInfo, winColors);
    
    if all(funcSuccess==1)
        success = 1;
        runOnceMore = false;
        break;
    elseif any(funcSuccess==-1)
        success = 0;
        runOnceMore = false;
        break
    end
            
    %if (GetSecs-loopTop)>params.waitForTolerance, warning('waitFor:tooSlow','waitFor exceeded latency tolerance - %s',datestr(now)); end; %warn tolerance exceeded -acs22dec2012
end



% for events that would have succeeded if things happen for *at least* a
% length of time, running the checker functions once more allows one to
% check if they have (since the loop will have broken out)
if runOnceMore
    loopTop = GetSecs;
    funcSuccess = zeros(1,length(eventCheckerFunction));
    for funcInd = 1:length(eventCheckerFunction)
        funcSuccess(funcInd) = eventCheckerFunction{funcInd}(loopStart, loopTop, eventFunctionInputs{funcInd}{:});
    end
    
    if all(funcSuccess)
        success = 1;
    elseif any(funcSuccess==-1)
        success = 0;
    end
end

% for funcInd = 1:length(eventCheckerFunction)
%     if isequal(eventCheckerFunction{funcInd}, @joystickAtariHold)
%         keyboard
%     end
% end

for funcInd = 1:length(eventCheckerFunction)
    clear(func2str(eventCheckerFunction{funcInd}));
end

