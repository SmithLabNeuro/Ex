function newStimPattern = closedLoopEpsGreedy(meanSpikeCount, ctlMsg, modelParams, expParams, binNum, controlCompSocket)

% persistent algo1 algo2 algo3 algo4
persistent algoStructs
% stimEffectPred errorEst exploredCnt stimEffectPredOrig exploredCntOrig errorEstOrig newStimChan currStimChan bciTrialCnt epsilon

ctlMsgChar = char(ctlMsg); 
% disp(ctlMsg)
% disp(ctlMsgChar)

% if ~isempty(ctlMsgChar) && ~strcmp(ctlMsgChar, 'ack') && ~(length(ctlMsg)==1) && (binNum==0)
if ~isempty(ctlMsgChar) && ~strcmp(ctlMsgChar, 'ack') && (binNum==0)
    % if (length(ctlMsg)==1) 
    while (length(ctlMsg)==1)||isempty(ctlMsg)
        disp('in')
        pause(0.01)
        ctlMsg = checkForControlMsgs(controlCompSocket);
        ctlMsg = uint8(ctlMsg);
    end
    
% if ~(isempty(ctlMsgChar)) && ~(strcmp(ctlMsgChar, 'ack')) && (binNum==0)
    ctlMsgDouble = typecast(ctlMsg, 'double');
    
    if isempty(ctlMsgDouble) % add this section to aviod an error ctlMsgDouble(1) exceeds the number of array elements
        updatedStimChan = 1;
        algoType = 1;
    else
        updatedStimChan = ctlMsgDouble(1) ;
        algoType = ctlMsgDouble(2); % 1=nostim, 2=random, 3=online only, 4=off+online
    end
    targetLatentDim = expParams.targetLatentDim;
    targetLatentVal = expParams.targetLatentVal;
    % targetLatentVal = mean(modelParams.posts) + 2*(std(modelParams.posts));
    stimChans = expParams.bciuStimChan;
    
    if isempty(algoStructs)
        stimEffectPredRand = modelParams.posts;
        stimEffectPredOff = modelParams.postsWithOff;
        exploredCnt = modelParams.exploredCntInCalib;
        exploredCntOff = modelParams.exploredCntInCalibOff;
        errorEst = abs(targetLatentVal - stimEffectPredRand);
        errorEstOff = abs(targetLatentVal - stimEffectPredOff);
        bciTrialCnt = 0;
        epsilon = expParams.epsilon;
        algo1 = struct('stimEffectPred', stimEffectPredRand, 'exploredCnt', exploredCnt, 'errorEst', errorEst, 'bciTrialCnt', bciTrialCnt, 'epsilon', epsilon);
        algo2 = struct('stimEffectPred', stimEffectPredRand, 'exploredCnt', exploredCnt, 'errorEst', errorEst, 'bciTrialCnt', bciTrialCnt, 'epsilon', epsilon);
        algo3 = struct('stimEffectPred', stimEffectPredRand, 'exploredCnt', exploredCnt, 'errorEst', errorEst, 'bciTrialCnt', bciTrialCnt, 'epsilon', epsilon);
        algo4 = struct('stimEffectPred', stimEffectPredOff, 'exploredCnt', exploredCntOff, 'errorEst', errorEstOff, 'bciTrialCnt', bciTrialCnt, 'epsilon', epsilon);
        algoStructs = struct('algo1', algo1, 'algo2', algo2, 'algo3', algo3, 'algo4', algo4);
    end
    

%     if isempty(stimEffectPred)
%         stimEffectPred = modelParams.posts;
%         exploredCnt = modelParams.exploredCntInCalib;
%         errorEst = abs(targetLatentVal - stimEffectPred);
%         bciTrialCnt = 0;
%         epsilon = expParams.epsilon;
%     end    

    % receive the stimulated  channel information
    % while ~matlabUDP2('check', controlCompSocket.receiver)
    %     pause(0.01)
    %     disp('wait')
    % end
    % receivedMsg = matlabUDP2('receive', controlCompSocket.receiver);
    % receivedMsgCasted = receivedMsg;
    % updatedStimChan = receivedMsgCasted(1);
    % save the original stim effect table
    
    if algoType==1
        algoStruct = algoStructs.algo1;
    elseif algoType==2
        algoStruct = algoStructs.algo2;
    elseif algoType==3
        algoStruct = algoStructs.algo3;
    elseif algoType==4
        algoStruct = algoStructs.algo4;
    end
    
    algoStruct.stimEffectPredOrig = algoStruct.stimEffectPred;
    algoStruct.exploredCntOrig = algoStruct.exploredCnt;
    algoStruct.errorEstOrig = algoStruct.errorEst;
    % targetLatentVal = expParams.targetLatentVal;

    if updatedStimChan~=0 && (algoType==3||algoType==4)
        % compute the mean posterior for the current trial data
        [zPost, ~, beta] = fastfa_estep(meanSpikeCount, modelParams.estParams);
        [zPostOrth, Lorth, TT] = orthogonalize(zPost.mean, modelParams.estParams.L);
        zPostMean = zPostOrth(targetLatentDim);

        if algoType==3
            % update the stim effect table
            if algoStruct.exploredCnt(updatedStimChan)==0
                algoStruct.stimEffectPred(updatedStimChan,:) = zPostMean;
            else
                algoStruct.stimEffectPred(updatedStimChan,:) = (algoStruct.stimEffectPred(updatedStimChan,:)*algoStruct.exploredCnt(updatedStimChan)+zPostMean')/(algoStruct.exploredCnt(updatedStimChan)+1);
            end
        elseif algoType==4
            lr = 0.1;
            currentPredError = zPostMean - algoStruct.stimEffectPred(updatedStimChan,:);
            algoStruct.stimEffectPred(updatedStimChan,:) = algoStruct.stimEffectPred(updatedStimChan) + lr*currentPredError;
        end
        % what is the objective here??
        % assuming one target value for now, and using absolute error
        algoStruct.errorEst(updatedStimChan,:) = abs(targetLatentVal - algoStruct.stimEffectPred(updatedStimChan,:)); 
        algoStruct.exploredCnt(updatedStimChan) = algoStruct.exploredCnt(updatedStimChan)+1;
    else
        zPostMean = 0;
    end
    
    algoStruct.bciTrialCnt = algoStruct.bciTrialCnt + 1;
    
    % decrease the value of epsilon for epsilon greedy algorithm
    if isfield(expParams, 'alpha')
        algoStruct.epsilon = algoStruct.epsilon * expParams.alpha;
        % disp(epsilon)
    end
    
    % pick up one uStim pattern
    % what is the algorithm of choosing one uStim pattern in greedy manner?
    if rand < algoStruct.epsilon
        indNewStimChan = randsample(stimChans,1);
        isGreedy = 0;
    else
        % indNewStimChan = find(errorEst==min(errorEst));
        indNewStimChan = find(sum(algoStruct.errorEst,2)==min(sum(algoStruct.errorEst,2)));
        if length(indNewStimChan)~=1
            indNewStimChan = randsample(indNewStimChan, 1);
        end
        isGreedy = 1;
    end
    % newStimChan = stimChans(indNewStimChan);
    newStimChan = indNewStimChan;
    
    % newStimPattern = [targetLatentVal,newStimChan];
    % newStimPattern = [newStimChan, isGreedy, targetLatentVal, zPostMean, bciTrialCnt];
    newStimPattern = [newStimChan, isGreedy, algoStruct.bciTrialCnt];
    
    if algoType==1
        algoStructs.algo1 = algoStruct;
    elseif algoType==2
        algoStructs.algo2 = algoStruct;
    elseif algoType==3
        algoStructs.algo3 = algoStruct;
    elseif algoType==4
        algoStructs.algo4 = algoStruct;
    end

    % save data
    subject = modelParams.subject;
    bciDecoderSaveDrive = '/home/smithlab/bciParameters/';
    bciDecoderRelativeSaveFolder = fullfile('closedLoopEpsGreedy/', subject, datestr(today, 'yymmdd'));
    bciDecoderSaveFolder = fullfile(bciDecoderSaveDrive, bciDecoderRelativeSaveFolder);
    success = mkdir(bciDecoderSaveFolder);
    if ~success
       fprintf('\nError creating new directory for BCI logs\n')
       fpringf('\nkeybord here...\n')
       keybord
    end

    % save the updated values to a file
    subjectCamelCase = lower(subject);
    subjectCamelCase(1) = upper(subjectCamelCase(1));
    bciDecoderSaveName = sprintf('%s%sclosedLoopEpsGreedy_%s', subjectCamelCase(1:2), datestr(today, 'yymmdd'), datestr(now, 'HH-MM-SS'));
    % save(fullfile(bciDecoderSaveFolder, bciDecoderSaveName), 'stimEffectPredOrig', 'exploredCntOrig', 'errorEstOrig', 'stimEffectPred', 'exploredCnt', 'errorEst', 'meanSpikeCount', 'newStimPattern', 'isGreedy', 'targetLatentVal', 'zPostMean', 'bciTrialCnt');
    save(fullfile(bciDecoderSaveFolder, bciDecoderSaveName), 'meanSpikeCount', 'newStimPattern', 'isGreedy', 'targetLatentVal', 'zPostMean', 'algoType', 'algoStructs', 'algoStruct');

    %sendMessageWaitAck(controlCompSocket, typecast(newStimChan, 'uint8'), 1000);
    % disp('output')
    disp(algoType)
    disp(newStimPattern)
    disp(sum(algoStruct.exploredCnt))
else
    newStimPattern = [];
end