function rtex(connectionType, inputBciParamFile)

global params codes
exGlobals;
xippmex('close');
disp('closed any existing xippmex connection')
%% setup xippmex connection for neural data
if nargin<1 || isempty(connectionType)
    connectionType = 'tcp';
end
if nargin<2
    inputBciParamFile = '';
end

status = false;
while ~status
    try
        status = xippmex(connectionType);
    catch error
        if contains(error.message, sprintf('xippmex command string ''%s'' not recognized', connectionType))
            status = xippmex;
            if ~status
                error('not sure what''s happening')
            end
        elseif contains(error.message, 'could not find NIP')
            xippmex('close')
        else
            rethrow(error)
        end
    end
    pause(0.2)
end

while ~status
    pause(0.2)
    status = xippmex(connectionType);
end

fprintf('xippmex running\n')


%% set up communication with control computer
controlCompSocket = initializeControlComputerCommunication();


%%
% controlCompSocket = 0;
microornano= 'nano';
okelecs = xippmex('elec',microornano);


while true
    % the below allows us to run a series of bciFunctions (say calibration
    % + increase control + final BCI) without having to rerun anything on
    % this computer
    [controlCompSocket, controlMsg] = checkIfControlComputerCommunicationNeedsReinitialization(controlCompSocket);

    if strcmp(controlMsg, 'readyBciParamFile')
        % get the bciParamFile; note that this loop is only meant to ping
        % the control computer for param files to run, once the function
        % associated with those param files ends, we assume a new one is
        % being expected.
        [bciParamFile, subject] = receiveBciDecoderInformation(controlCompSocket);
    else
        bciParamFile = '';
        subject = '';
    end
    
    if ~isempty(inputBciParamFile)
        if ~isempty(bciParamFile)
            warning('using param file sent by control computer in lieu of input in BCI computer');
        else
            bciParamFile = inputBciParamFile;
        end
    end
    if ~isempty(bciParamFile)
        if ~exist(bciParamFile, 'file')
            decoderParameterLocation = params.bciDecoderBasePathBciComputer;
            xmlParameterFolder = params.bciDecoderXmlParamFolder;
            bciParamFile = fullfile(decoderParameterLocation, xmlParameterFolder, bciParamFile);
        end
        [~,machineInit] = system('hostname');
        machine = lower(deblank(cell2mat(regexp(machineInit,'^[^\.]+','match'))));
        [~, expParams, ~, ~] = readExperiment(bciParamFile, subject, machine);
        % Save the subject in the struct
        expParams.subject = subject;
        bciWrapper = str2func([expParams(1).bciStyle 'Bci']);
        % bciFunction does ALL the work, and when it finished it pops us
        % back out here
%         bciFunction = str2func(expParams(1).exFileName);
        bciWrapper(controlCompSocket, expParams, okelecs)
        disp('BCI function done running, waiting for next one')
    end
    
    pause(1);
end

end



