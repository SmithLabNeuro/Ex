function bciEnd = checkIfBciEnd(controlCompSocket)
    % checks if 'bciEnd' is received to tell us if a new decoder is needed
    controlMsg = checkForControlMsgs(controlCompSocket);
    if strcmp(controlMsg, 'bciEnd')
        bciEnd = true;
    else
        bciEnd = false;
    end
end