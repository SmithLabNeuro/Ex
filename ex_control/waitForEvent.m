function [success, funcSuccess] = waitForEvent(timeout, eventCheckerFunction, eventFunctionInputs)
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
while (toc(ticStart)*1000) <= timeout
    loopTop = GetSecs;
    
    funcSuccess = zeros(1,length(eventCheckerFunction));
    for funcInd = 1:length(eventCheckerFunction)
        funcSuccess(funcInd) = eventCheckerFunction{funcInd}(loopStart, loopTop, eventFunctionInputs{funcInd}{:});
    end
    
    if all(funcSuccess==1)
        success = 1;
        runOnceMore = false;
        break;
    elseif any(funcSuccess==-1)
        success = 0;
        runOnceMore = false;
        break
    end
    
    % let's one press e.g. 'q' in order to quit a stimulus in the middle
    % with runex
    if keyboardEvents()
        success = 0;
        runOnceMore = false;
        break;
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

