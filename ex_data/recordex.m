function recordex(connectionType)

global params codes
% do this before adding paths
exGlobals; % grab the global parameters

% addpath(genpath('C:\Users\smithlab\Ex\ex_control'));
% addpath(genpath('C:\Users\smithlab\spikesort\'))

if nargin<1
    connectionType = 'tcp';
end


% The first close (at least under UDP), seems to recognize a connection
% upon Matlab start (that existed from before Matlab started!) but not
% close it--this recognition makes the attempt at opening fail, because
% xippmex is in the wrong state to open new connections. As a result,
% recordex fails. In any case, a second 'close' seems to fix this...
% xippmex('close'); BUT! Things get more confusing, so see the comment in
% the try-catch below.
try
    % In a twist of confusion, when you assign an output to
    % xippmex('close'), xippmex actually errors if everything is fine
    % because apparently it doesn't assign an output when things are
    % correctly closed. If it *does* assign an output, apparently that's
    % an indication that things didn't close correctly, at which point we
    % need to close again. So, here, if an output is assigned it continues
    % to the second xippmex('close'), and if an output is not assigned but
    % we request it, it errors and as long as that output error was what
    % happened, recordex continues happily.
    statCheck = xippmex('close');
    pause(0.1); % pause needed I think to prevent crash
    xippmex('close')
catch error
    if ~contains(error.message, 'One or more output arguments not assigned during call to "xippmex".')
        rethrow(error)
    end
end

% go = waitForMessage;
status = false;
while ~status
    pause(0.2)
    try
        status = xippmex(connectionType);
    catch error
        if contains(error.message, sprintf('xippmex command string ''%s'' not recognized', connectionType))
            status = xippmex(connectionType);
            if ~status
                error('not sure what''s happening')
            end
        elseif ~contains(error.message, 'could not find NIP')
            rethrow(error)
        end
    end
end

% have to add an operator if this is a TCP connection
if strcmp(connectionType, 'tcp')
    operatorId = 129; % this it the last octet of the IP address of Ripple
    timeInLoop = 1;
    while true
        try
            xippmex('addoper', operatorId);
        catch err
            if isempty(strfind(err.message, 'operator not found'))
                rethrow(err)
            else
                % for some reason if it's running right as Trellis opens
                % this error appears, so we make sure to only show it if it
                % continually appears...
                if timeInLoop > 1 
                    fprintf('Trouble connecting to Trellis. Is it open?\n')
                end
                pause(0.2)
                timeInLoop = timeInLoop+1;
                continue % this keeps us in the loop
            end
        end
        % don't like this construction but it pops us out of the loop if
        % the connection succeeds
        fprintf('Connected to Trellis.\n')
        break 
    end
end

fprintf('xippmex running\n')

% open communication port to control computer; numbers are selected based
% on ports in exGlobals; at some point would be ideal to move these values
% into there
udpPortFromControlReceive = 4245;
udpPortFromControlSend = 4246;
controlIPAddress = '128.2.245.223';
maxNumBytes = 2000; % this value mimics matlabUDP2.h; should be enough
udpr = dsp.UDPReceiver('RemoteIPAddress', controlIPAddress, 'LocalIPPort', udpPortFromControlReceive, 'MaximumMessageLength', maxNumBytes);
udps = dsp.UDPSender('RemoteIPAddress', controlIPAddress, 'RemoteIPPort', udpPortFromControlSend);

socketsControlComm.receiver = udpr;
socketsControlComm.sender = udps;

defaultPauseTime = .1;
longPauseTime = 1;
pauseTime = defaultPauseTime;
dataPath = 'E:\';
bciDecodersBasePath = params.bciDecoderBasePathDataComputer;
bciDecodersParameterPath = fullfile(bciDecodersBasePath, params.bciDecoderXmlParamFolder);
while true
    
    msg = receiveMessageSendAck(socketsControlComm);
%     if isempty(msg)
%         continue
    if strcmp(msg, 'record')
        subjectName = receiveMessageSendAck(socketsControlComm);
        sessionNum = receiveMessageSendAck(socketsControlComm);
        xmlName = receiveMessageSendAck(socketsControlComm);
        
        % I call it base because Trellis autoincrememnts...
        filenameBase = [upper(subjectName(1)) lower(subjectName(2)) datestr(today(), 'YYmmdd') '_s' sessionNum '_' xmlName '_'];
        fileLoc = fullfile(dataPath, lower(subjectName), filenameBase);
        if isempty(dir([fileLoc, '*']))
            incSt = 1;
        else
            fileInfoExist = dir([fileLoc, '*']);
            fileNames = {fileInfoExist.name};
            fileNameNumLocs = cellfun(@(fn) regexp(fn, '\d{4}'), fileNames, 'uni', 0);
            % this will also capture the date, so we grab the last index
            % which is the increment number
            fileNameIncNums = cellfun(@(fn, numLoc) str2double(fn(numLoc(end)+(0:3))), fileNames, fileNameNumLocs);
            incSt = max(fileNameIncNums)+1;
        end
        % the 0 says don't stop, the 1 says autoincrement, and the incSt
        % tells it the current increment (or if empty keeps whatever's
        % there)
        disp(fileLoc)
        refFlPath = sprintf('%s%04d', fileLoc(length(dataPath)+1:end), incSt);
        sendMessageWaitAck(socketsControlComm, uint8(refFlPath));
        recordingInfo = xippmex('trial','recording',fileLoc, 0, 1, incSt);
        if ~strcmp(recordingInfo.status, 'recording')
            pause(0.1)
            recordingInfo = xippmex('trial','recording',fileLoc, 0, 1, incSt);
            if ~strcmp(recordingInfo.status, 'recording')
                keyboard
                sendMessageWaitAck(socketsControlComm, uint8(recordingInfo.status));
            else
                sendMessageWaitAck(socketsControlComm, uint8(recordingInfo.status));
            end
        else
            sendMessageWaitAck(socketsControlComm, uint8(recordingInfo.status));
        end
        pauseTime = longPauseTime;
    elseif strcmp(msg, 'stopRecording')
        fileSavedInfo = xippmex('trial','stopped');
        relFlPath = [fileSavedInfo.filebase(length(dataPath)+1:end) sprintf('%04d', fileSavedInfo.incr_num)];
%         sendMessageWaitAck(socketsControlComm, uint8(relFlPath));
%         copyToRaptorRigM;
%         copyToRaptor;
%         disp('Copying done')
        pauseTime = defaultPauseTime;
    elseif strcmp(msg, 'pauseRecording')
        fileSavedInfo = xippmex('trial','paused')
        relFlPath = [fileSavedInfo.filebase(length(dataPath)+1:end) sprintf('%04d', fileSavedInfo.incr_num-1)];
        sendMessageWaitAck(socketsControlComm, uint8(relFlPath));
%         copyToRaptorRigM;
%         copyToRaptor;
%         disp('Copying done')
        pauseTime = defaultPauseTime;
    elseif strcmp(msg, 'trainDecoder') || strcmp(msg, 'trainDecoderPause')
        % the code below is for training a BCI decoder when requested by
        % the control computer.
        % Stops recording for trellis
        if strcmp(msg, 'trainDecoder')
            fileSavedInfo = xippmex('trial','stopped')
            % Checks that recording stopped and if so, keeps track of increment
            % suffix
            if strcmp(fileSavedInfo.status, 'stopped')
                incSaved = fileSavedInfo.incr_num;
                % Gets the fileName of nev file that was just made
                relFlPath = [fileSavedInfo.filebase(length(dataPath)+1:end) sprintf('%04d', incSaved)];
                if isempty(dir([fileSavedInfo.filebase sprintf('%04d', incSaved) '.nev']))
                    disp('Trellis decided to fake us out on the increment number...')
                    relFlPath = [fileSavedInfo.filebase(length(dataPath)+1:end) sprintf('%04d', incSaved-1)];
                end
            else
                fprintf(['\nIMPORTANT: You are in debug mode.\n' ...
                    'Try running xippmex(''trial'',''stopped'' command above and then press F5 if status returns ''stopped''\n']);
                keyboard
                incSaved = fileSavedInfo.incr_num;
            end
        elseif strcmp(msg, 'trainDecoderPause')
            fileSavedInfo = xippmex('trial','paused')
        end
        % Receive filename to use to train decoder from control computer
        % (decoderCalibrationFunctionName)
        decoderTrainFunctionName = receiveMessageSendAck(socketsControlComm);
        decoderTrainFunction = str2func(decoderTrainFunctionName);
        % Receives filename for parameters (file is in a drive shared with
        % between bci/data computer to make it easier to ensure that
        % decoder training and decoder application parameters are the same)
        % to train decoder via bciDecoderParamFile from control
        decoderTrainParameterFile = receiveMessageSendAck(socketsControlComm);
        decoderTrainParameterFilepath = fullfile(bciDecodersParameterPath, decoderTrainParameterFile);
        % Reads in xml files for decoder parameters to cell structs
        [~,machineInit] = system('hostname');
        machine = lower(deblank(cell2mat(regexp(machineInit, '^[^\.]+', 'match'))));
        [~, trainParams, ~, ~] = readExperiment(decoderTrainParameterFilepath, '',machine);
        
        % send over the BCI parameters to Control computer so they can be
        % saved to the NEV
        sendMessageWaitAck(socketsControlComm, 'startSendingAsciiParameters');
        sendStructAsAscii(trainParams, socketsControlComm);
        sendMessageWaitAck(socketsControlComm, 'endSendingAsciiParameters');
        % Receives nev file names to calibrate on from control computer
        nevBaseRelFilepath = receiveMessageSendAck(socketsControlComm);
        nevBaseFullFilepath = fullfile(dataPath, nevBaseRelFilepath);
        nevFilesForTraining = receiveMessageSendAck(socketsControlComm);
        sepTrainRelFilepaths = strsplit(nevFilesForTraining, '\n');
        sepTrainFullFilepaths = cellfun(@(rfp) fullfile(dataPath, rfp), sepTrainRelFilepaths, 'uni', 0);

        % Train decoder on the files received from the control computer
        trainedDecoderInfoChar = decoderTrainFunction(socketsControlComm, nevBaseFullFilepath, sepTrainFullFilepaths, trainParams, subjectName);
        % Send trained decoder info to the control computer as a char (i.e.
        % the location of the decoder, or the location of two decoders
        % split by a newline)
        sendMessageWaitAck(socketsControlComm, uint8(trainedDecoderInfoChar));
        
        if strcmp(msg, 'trainDecoderPause')
            recordingInfo = xippmex('trial','paused') % this is how you restart a pause?
            if ~strcmp(recordingInfo.status, 'recording')
                pause(1)
                recordingInfo = xippmex('trial','paused')
                if ~strcmp(recordingInfo.status, 'recording')
                    keyboard
                    %last try?
                    recordingInfo = xippmex('trial','paused')
                    sendMessageWaitAck(socketsControlComm, uint8('recording'));
                else
                    sendMessageWaitAck(socketsControlComm, uint8(recordingInfo.status));
                end
            else
                sendMessageWaitAck(socketsControlComm, uint8(recordingInfo.status));
            end
        end
        
        pauseTime = longPauseTime;
    elseif strcmp(msg, 'sendDecoderParameters')
        filePausedInfo = xippmex('trial','paused')
        
        decoderTrainParameterFile = receiveMessageSendAck(socketsControlComm);
        decoderTrainParameterFilepath = fullfile(bciDecodersParameterPath, decoderTrainParameterFile);
        [~,machineInit] = system('hostname');
        machine = lower(deblank(cell2mat(regexp(machineInit, '^[^\.]+', 'match'))));
        [~, trainParams, ~, ~] = readExperiment(decoderTrainParameterFilepath, '',machine);
        
        % send over the BCI parameters so they can be saved to the NEV
        sendMessageWaitAck(socketsControlComm, 'startSendingAsciiParameters');
        sendStructAsAscii(trainParams, socketsControlComm);
        sendMessageWaitAck(socketsControlComm, 'endSendingAsciiParameters');
        
        recordingInfo = xippmex('trial','paused') % this is how you restart a pause?
        if ~strcmp(recordingInfo.status, 'recording')
            pause(1)
            recordingInfo = xippmex('trial','paused')
            if ~strcmp(recordingInfo.status, 'recording')
                keyboard
                %last try?
                recordingInfo = xippmex('trial','paused')
                sendMessageWaitAck(socketsControlComm, uint8('recording'));
            else
                sendMessageWaitAck(socketsControlComm, uint8(recordingInfo.status));
            end
        else
            sendMessageWaitAck(socketsControlComm, uint8(recordingInfo.status));
        end
    elseif strcmp(msg, 'restartPausedRecording')
        recordingInfo = xippmex('trial','paused') % this is how you restart a pause?
        if ~strcmp(recordingInfo.status, 'recording')
            pause(1)
            recordingInfo = xippmex('trial','paused')
            if ~strcmp(recordingInfo.status, 'recording')
                keyboard
                %last try?
                recordingInfo = xippmex('trial','paused')
                sendMessageWaitAck(socketsControlComm, uint8('recording'));
            else
                sendMessageWaitAck(socketsControlComm, uint8(recordingInfo.status));
            end
        else
            sendMessageWaitAck(socketsControlComm, uint8(recordingInfo.status));
        end
    elseif strcmp(msg, 'sessionEnd')
        addpath 'C:\Users\rigmdata\Documents\labcode\rigutils'
        try
            copyToRaptorRigM;
        catch
            copyToRaptor;
        end
        close all;
%         break;
    end
    
    pause(pauseTime);
end
end
