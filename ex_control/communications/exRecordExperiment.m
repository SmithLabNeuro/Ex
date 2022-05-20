%% function to initialize recording on the data computer using appropriate parameters
    function [isRecording, defaultRunexPrompt] = exRecordExperiment(socketsDatComp, isRecording, sessionInfo, xmlParams, outfilename, defaultRunexPrompt)
    global trialData notes
        timeout = 10;
        
        promptSt = 'Communicating with data computer to start recording...';
        trialData{4} = promptSt;
        drawTrialData();
        
%         if false;
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
%             neuralOutName = receiveMessageSendAck(socketsDatComp, timeout);
%             sqlDb.insert('experiment_info', {'animal', 'date', 'task', 'rig', 'behavior_output_name', 'neural_output_name'}, {params.SubjectID, datestr(today), xmlBase, params.machine, outfilename, neuralOutName})
            currPrompt = 'stop (r)ecording neural data';
            nextPrompt = '(r)ecord neural data';
        end
%         end
        
        promptSt = strfind(defaultRunexPrompt, currPrompt);
        defaultRunexPrompt(promptSt:promptSt+length(currPrompt)-1) = [];
        defaultRunexPrompt = [defaultRunexPrompt(1:promptSt-1), nextPrompt, defaultRunexPrompt(promptSt:end)];
        isRecording = ~isRecording;
        trialData{4} = defaultRunexPrompt;
        drawTrialData();

    end
    