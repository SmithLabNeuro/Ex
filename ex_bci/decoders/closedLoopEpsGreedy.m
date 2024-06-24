function newStimPattern = closedLoopEpsGreedy(meanSpikeCount, ctlMsg, modelParams, expParams, binNum, controlCompSocket)

% persistent algo1 algo2 algo3 algo4
persistent algoStructs bciTrialCntOverall
% stimEffectPred errorEst exploredCnt stimEffectPredOrig exploredCntOrig errorEstOrig newStimChan currStimChan bciTrialCnt epsilon

ctlMsgChar = char(ctlMsg); 
% disp(ctlMsg)
% disp(ctlMsgChar)

% if ~isempty(ctlMsgChar) && ~strcmp(ctlMsgChar, 'ack') && ~(length(ctlMsg)==1) && (binNum==0)
in_counter = 0;
if ~isempty(ctlMsgChar) && ~strcmp(ctlMsgChar, 'ack') && (binNum==0)
    % if (length(ctlMsg)==1) 
    while (length(ctlMsg)==1)||isempty(ctlMsg)
        disp('in')
        pause(0.01)
        ctlMsg = checkForControlMsgs(controlCompSocket);
        ctlMsg = uint8(ctlMsg);
        in_counter = in_counter + 1;
        if in_counter > 100
            break
        end
    end
    
% if ~(isempty(ctlMsgChar)) && ~(strcmp(ctlMsgChar, 'ack')) && (binNum==0)
    if in_counter < 100
        ctlMsgDouble = typecast(ctlMsg, 'double');

        updatedStimChan = ctlMsgDouble(1) ;
        algoType = ctlMsgDouble(2); % 1=nostim, 2=random, 3=online only, 4=off+online
    else % if observed too many in
        updatedStimChan = 0;
        algoType = 1; % 1=nostim, 2=random, 3=online only, 4=off+online
    end
%     if isempty(ctlMsgDouble) % add this section to aviod an error ctlMsgDouble(1) exceeds the number of array elements
%         updatedStimChan = 1;
%         algoType = 1;
%     else
%         updatedStimChan = ctlMsgDouble(1) ;
%         algoType = ctlMsgDouble(2); % 1=nostim, 2=random, 3=online only, 4=off+online
%     end

    % targetLatentVal = mean(modelParams.posts) + 2*(std(modelParams.posts));
    targetLatentDim = expParams.targetLatentDim;
    
    if isfield(expParams, 'bciTargetSwitchTrial')
        if bciTrialCntOverall>=expParams.bciTargetSwitchTrial
            targetLatentVal = [0.3, 0];
        else
            targetLatentVal = expParams.targetLatentVal;
        end
    else
        targetLatentVal = expParams.targetLatentVal;
    end
    
    if isempty(algoStructs)
        stimEffectPredRand = modelParams.posts;
        stimEffectPredOff = modelParams.postsWithOff;
        exploredCnt = modelParams.exploredCntInCalib;
        exploredCntOff = modelParams.exploredCntInCalibOff;
        % errorEst = abs(targetLatentVal - stimEffectPredRand);
        % errorEst = mean((targetLatentVal - stimEffectPredRand).^2, 2); %L2
        errorEst = mean(abs(targetLatentVal - stimEffectPredRand), 2); %L1
        % errorEst = (abs(targetLatentVal(1) - stimEffectPredRand(:,1))+0.5*abs(targetLatentVal(2) - stimEffectPredRand(:,2)))/2; %L1 with weight
        % errorEstOff = abs(targetLatentVal - stimEffectPredOff);
        % errorEstOff = mean((targetLatentVal - stimEffectPredOff).^2, 2);
        errorEstOff = mean(abs(targetLatentVal - stimEffectPredOff), 2); %L1
        % errorEstOff = (abs(targetLatentVal(1) - stimEffectPredOff(:,1))+0.5*abs(targetLatentVal(2) - stimEffectPredOff(:,2)))/2; %L1 with weight
        bciTrialCnt = 0;
        bciTrialCntOverall = 0;
        epsilon = expParams.epsilon;
        algo1 = struct('stimEffectPred', stimEffectPredRand, 'exploredCnt', exploredCnt, 'errorEst', errorEst, 'bciTrialCnt', bciTrialCnt, 'epsilon', epsilon);
        algo2 = struct('stimEffectPred', stimEffectPredRand, 'exploredCnt', exploredCnt, 'errorEst', errorEst, 'bciTrialCnt', bciTrialCnt, 'epsilon', epsilon);
        algo3 = struct('stimEffectPred', stimEffectPredRand, 'exploredCnt', exploredCnt, 'errorEst', errorEst, 'bciTrialCnt', bciTrialCnt, 'epsilon', epsilon); % single, model-free
        % algo3 = struct('stimEffectPred', stimEffectPredOff, 'exploredCnt', exploredCntOff, 'errorEst', errorEstOff, 'bciTrialCnt', bciTrialCnt, 'epsilon', epsilon); % double, model-base 
        algo4 = struct('stimEffectPred', stimEffectPredOff, 'exploredCnt', exploredCntOff, 'errorEst', errorEstOff, 'bciTrialCnt', bciTrialCnt, 'epsilon', epsilon);
        algoStructs = struct('algo1', algo1, 'algo2', algo2, 'algo3', algo3, 'algo4', algo4);
    end
    
    if isfield(expParams, 'bciTargetSwitchTrial')
        if bciTrialCntOverall==expParams.bciTargetSwitchTrial
            algoStructs.algo3.errorEst = mean(abs(targetLatentVal - algoStructs.algo3.stimEffectPred), 2); %L1
            % algoStructs.algo3.errorEst = (abs(targetLatentVal(1) - algoStructs.algo3.stimEffectPred(:,1))+0.5*abs(targetLatentVal(2) - algoStructs.algo3.stimEffectPred(:,2)))/2; %L1 with weight
            algoStructs.algo4.errorEst = mean(abs(targetLatentVal - algoStructs.algo4.stimEffectPred), 2); %L1
            % algoStructs.algo4.errorEst = (abs(targetLatentVal(1) - algoStructs.algo4.stimEffectPred(:,1))+0.5*abs(targetLatentVal(2) - algoStructs.algo4.stimEffectPred(:,2)))/2; %L1 with weight
        end
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
    valid_index_range = 1:length(algoStruct.stimEffectPred);

    if updatedStimChan~=0 && (algoType==3||algoType==4)
        % compute the mean posterior for the current trial data
        [zPost, ~, ~] = fastfa_estep(meanSpikeCount, modelParams.estParams);   
        zPostOrth = zPost.mean;
        
        %[zPost, ~, ~] = fastfa_estep(meanSpikeCount, modelParams.estParams);
        %[zPostOrth, Lorth, ~] = orthogonalize(zPost.mean, modelParams.estParams.L);
        zPostMean = zPostOrth(targetLatentDim);

        if algoType==3
            % update the stim effect table
            if algoStruct.exploredCnt(updatedStimChan)==0
                algoStruct.stimEffectPred(updatedStimChan,:) = zPostMean;
            else
                % v2) update with average
                % algoStruct.stimEffectPred(updatedStimChan,:) = (algoStruct.stimEffectPred(updatedStimChan,:)*algoStruct.exploredCnt(updatedStimChan)+zPostMean')/(algoStruct.exploredCnt(updatedStimChan)+1);
                % v1) udate with learning rate
                lr_algo3_clip = 0.1;
                lr_algo3 = max(lr_algo3_clip, 1/(algoStruct.exploredCnt(updatedStimChan)+1));
                currentPredError = zPostMean - algoStruct.stimEffectPred(updatedStimChan,:)';
                algoStruct.stimEffectPred(updatedStimChan,:) = algoStruct.stimEffectPred(updatedStimChan,:) + lr_algo3*currentPredError';            
            end
            if isfield(expParams, 'filterValidChan')
                if expParams.filterValidChan
                    valid_index_range = 1:96;
                end
            end
            
        elseif algoType==4
            % v1) update with learning rate
            lr_algo4_clip = 0.1;
            lr_algo4 = max(lr_algo4_clip, 1/(algoStruct.exploredCnt(updatedStimChan)+2));
            currentPredError = zPostMean - algoStruct.stimEffectPred(updatedStimChan,:)';
            algoStruct.stimEffectPred(updatedStimChan,:) = algoStruct.stimEffectPred(updatedStimChan,:) + lr_algo4*currentPredError';
            if isfield(expParams, 'filterValidChan')
                if expParams.filterValidChan
                    valid_index_range = 97:length(algoStruct.stimEffectPred);
                end
            end
            % v2) update with average
            % algoStruct.stimEffectPred(updatedStimChan,:) = (algoStruct.stimEffectPred(updatedStimChan,:)*algoStruct.exploredCnt(updatedStimChan)+zPostMean')/(algoStruct.exploredCnt(updatedStimChan)+1);
            % use below when exploreation Cnt in trial Type 4 after
            % calibration is 0
            % algoStruct.stimEffectPred(updatedStimChan,:) = (algoStruct.stimEffectPred(updatedStimChan,:)*(algoStruct.exploredCnt(updatedStimChan)+1)+zPostMean')/(algoStruct.exploredCnt(updatedStimChan)+2);
        end
        % what is the objective here??
        % assuming one target value for now, and using absolute error
        % algoStruct.errorEst(updatedStimChan,:) = abs(targetLatentVal - algoStruct.stimEffectPred(updatedStimChan,:));
        % multi dimensional target value with mean squared error
        % algoStruct.errorEst(updatedStimChan,:) = mean((targetLatentVal - algoStruct.stimEffectPred(updatedStimChan,:)).^2, 2); % L2
        algoStruct.errorEst(updatedStimChan,:) = mean(abs(targetLatentVal - algoStruct.stimEffectPred(updatedStimChan,:)), 2); % L1
        % algoStruct.errorEst(updatedStimChan,:) = (abs(targetLatentVal(1) - algoStruct.stimEffectPred(updatedStimChan,1))+0.5*abs(targetLatentVal(2) - algoStruct.stimEffectPred(updatedStimChan,2)))/2; %L1 with weight
        algoStruct.exploredCnt(updatedStimChan) = algoStruct.exploredCnt(updatedStimChan)+1;
    else
        zPostMean = 0;
    end
    
    algoStruct.bciTrialCnt = algoStruct.bciTrialCnt + 1;
    bciTrialCntOverall = bciTrialCntOverall + 1;
    
    % decrease the value of epsilon for epsilon greedy algorithm
    if isfield(expParams, 'alpha')
        algoStruct.epsilon = algoStruct.epsilon * expParams.alpha;
        % disp(epsilon)
    end
    
    % pick up one uStim pattern
    % what is the algorithm of choosing one uStim pattern in greedy manner?
    if rand < algoStruct.epsilon
        indNewStimChan = randsample(valid_index_range,1); % choose a random index corresponding to each uStim pattern
        isGreedy = 0;
    else
        % indNewStimChan = find(errorEst==min(errorEst));
        indNewStimChan = find(sum(algoStruct.errorEst,2)==min(sum(algoStruct.errorEst(valid_index_range,:),2)));
        if length(indNewStimChan)~=1
            indNewStimChan = randsample(indNewStimChan, 1);
        end
        isGreedy = 1;
    end

    if in_counter > 100
        isGreedy = 2; % something wrong in this trial so don't want to update nextStimParam on control
    end
    newStimChan = indNewStimChan;
    
    % newStimPattern = [targetLatentVal,newStimChan];
    % newStimPattern = [newStimChan, isGreedy, targetLatentVal, zPostMean, bciTrialCnt];
    % newStimPattern = [newStimChan, isGreedy, algoStruct.bciTrialCnt];
    newStimPattern = [newStimChan, isGreedy, bciTrialCntOverall];
    
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
    save(fullfile(bciDecoderSaveFolder, bciDecoderSaveName), 'meanSpikeCount', 'newStimPattern', 'isGreedy', 'targetLatentVal', 'zPostMean', 'algoType', 'algoStructs', 'algoStruct', 'bciTrialCntOverall');

    %sendMessageWaitAck(controlCompSocket, typecast(newStimChan, 'uint8'), 1000);
    % disp('output')
    disp(algoType)
    disp(newStimPattern)
    disp(sum(algoStruct.exploredCnt))
else
    newStimPattern = [];
end