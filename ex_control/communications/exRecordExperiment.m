%% function to initialize recording on the data computer using appropriate parameters
    function [isRecording, defaultRunexPrompt] = exRecordExperiment(socketsDatComp, isRecording, xmlParams, outfilename, defaultRunexPrompt)
    global trialData wins
        timeout = 10;
        
        promptSt = 'Communicating with data computer to start recording...';
        trialData{wins.trialData.errorLine} = promptSt;
        drawTrialData();
        
        if ~isRecording
            xmlFile = xmlParams.xmlFile;
            recStarted = startNevRecording(socketsDatComp, xmlFile, outfilename);
            
            if recStarted
                currPrompt = '(r)ecord neural data';
                nextPrompt = 'stop (r)ecording neural data';
            else
                return
            end
        else
            rc = sendMessageWaitAck(socketsDatComp, 'stopRecording', timeout);
            currPrompt = 'stop (r)ecording neural data';
            nextPrompt = '(r)ecord neural data';
        end
        
        promptSt = strfind(defaultRunexPrompt, currPrompt);
        defaultRunexPrompt(promptSt:promptSt+length(currPrompt)-1) = [];
        defaultRunexPrompt = [defaultRunexPrompt(1:promptSt-1), nextPrompt, defaultRunexPrompt(promptSt:end)];
        isRecording = ~isRecording;
        trialData{wins.trialData.errorLine} = defaultRunexPrompt;
        drawTrialData();

    end
    