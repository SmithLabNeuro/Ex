function recordex(connectionType)

global params codes
% do this before adding paths
exGlobals; % grab the global parameters

addpath(genpath('C:\Users\rigmdata\Documents\Ex\ex_control'));
addpath(genpath('C:\Users\rigmdata\spikesort\'))

xippmex('close');
if nargin<1
    connectionType = 'tcp';
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

pauseTime = 0.1;
dataPath = 'E:\';
bciDecodersBasePath = params.bciDecoderBasePathDataComputer;
bciDecodersParameterPath = fullfile(bciDecodersBasePath, params.bciDecoderXmlParamFolder);
while true
    
    msg = waitForMessage(udpr, udps);
%     if isempty(msg)
%         continue
    if strcmp(msg, 'record')
        subjectName = waitForMessage(udpr, udps);
        sessionNum = waitForMessage(udpr, udps);
        xmlName = waitForMessage(udpr, udps);
        
        % I call it base because Trellis autoincrememnts...
        filenameBase = [upper(subjectName(1)) lower(subjectName(2)) datestr(today(), 'YYmmdd') '_s' sessionNum 'a_' xmlName '_'];
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
        pauseTime = 1;
    elseif strcmp(msg, 'stopRecording')
        fileSavedInfo = xippmex('trial','stopped');
        relFlPath = [fileSavedInfo.filebase(length(dataPath)+1:end) sprintf('%04d', fileSavedInfo.incr_num)];
%         sendMessageWaitAck(socketsControlComm, uint8(relFlPath));
%         copyToRaptorRigM;
%         copyToRaptor;
%         disp('Copying done')
        pauseTime = 0.1;
    elseif strcmp(msg, 'pauseRecording')
        fileSavedInfo = xippmex('trial','paused')
        relFlPath = [fileSavedInfo.filebase(length(dataPath)+1:end) sprintf('%04d', fileSavedInfo.incr_num-1)];
        sendMessageWaitAck(socketsControlComm, uint8(relFlPath));
%         copyToRaptorRigM;
%         copyToRaptor;
%         disp('Copying done')
        pauseTime = 0.1;
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
        decoderTrainFunctionName = waitForMessage(udpr, udps);
        decoderTrainFunction = str2func(decoderTrainFunctionName);
        % Receives filename for parameters (file is in a drive shared with
        % between bci/data computer to make it easier to ensure that
        % decoder training and decoder application parameters are the same)
        % to train decoder via bciDecoderParamFile from control
        decoderTrainParameterFile = waitForMessage(udpr, udps);
        decoderTrainParameterFilepath = fullfile(bciDecodersParameterPath, decoderTrainParameterFile);
        % Reads in xml files for decoder parameters to cell structs
        [~, trainParams, ~, ~] = readExperiment(decoderTrainParameterFilepath, '');
        
        % send over the BCI parameters to Control computer so they can be
        % saved to the NEV
        sendMessageWaitAck(socketsControlComm, 'startSendingAsciiParameters');
        sendStructAsAscii(trainParams, socketsControlComm);
        sendMessageWaitAck(socketsControlComm, 'endSendingAsciiParameters');
        % Receives nev file names to calibrate on from control computer
        nevBaseRelFilepath = waitForMessage(udpr, udps);
        nevBaseFullFilepath = fullfile(dataPath, nevBaseRelFilepath);
        nevFilesForTraining = waitForMessage(udpr, udps);
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
        
        pauseTime = 1;
    elseif strcmp(msg, 'sendDecoderParameters')
        filePausedInfo = xippmex('trial','paused')
        
        decoderTrainParameterFile = waitForMessage(udpr, udps);
        decoderTrainParameterFilepath = fullfile(bciDecodersParameterPath, decoderTrainParameterFile);
        [~, trainParams, ~, ~] = readExperiment(decoderTrainParameterFilepath, '');
        
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
        copyToRaptorRigM;
        close all;
%         break;
    end
    
    pause(pauseTime);
end
end

% function msg = waitForMessage(udpr, udps)
%     msg = char(udpr())';
%     while isempty(msg)
%         pause(0.1)
%         msg = char(udpr())';
%         if length(msg) == udpr.MaximumMessageLength
%             warning('a message longer than the maximum buffer length may have been received; if so, this can be changed for udpr with the MaximumMessageLength parameter so there''s hope');
%         end
%     end
% %     disp(msg)
%     udps(uint8('ack'));
% end
% 
% function msg = sendMessageWaitAck(udps, udpr, msgToSend)
%     if ischar(msgToSend)
%         msgToSend = uint8(msgToSend);
%     end
% 
%     udps(msgToSend);
%     msg = char(udpr());
%     while isempty(msg)
%         pause(0.1)
%         msg = char(udpr())';
%     end
%     
%     if ~strcmp(msg, 'ack')
%         keyboard
%     end
% 
% end
% 
% function sendStructAsAscii(trainParams, udpr, udps)
% 
% %function sendStructAsAscii(s)
% %
% % Send a struct as ascii over udp
% %
% 
% fields = fieldnames(trainParams);
% 
% for i = 1:length(fields)
%     val = trainParams.(fields{i});
%     if iscell(val) % works for single layer cells, not nested ones
%         for clInd = 1:length(val)
%             valC = num2str(val{clInd}(:)'); % needs to be a row for the line below
%             m = [fields{i} '{' num2str(clInd) '}=' valC ';'];
%             
%             sendMessageWaitAck(udps, udpr, m);
%         end
%     else
%         val = num2str(val(:)'); % needs to be a row for the line below
%         m = [fields{i} '=' val ';'];
%         
%         sendMessageWaitAck(udps, udpr, m);
%     end
%     
% 
% end
% 
% end