function boundedBci(controlCompSocket, expParams, okelecs)

global params codes
if isfield(expParams, "saveOnlineMatFile")
    saveOnlineMatFile = expParams.saveOnlineMatFile; 
else
    saveOnlineMatFile = false;
end
if isfield(expParams, "refreshOutputEachTrial")
    refreshOutput = expParams.refreshOutputEachTrial;
else
    refreshOutput = true; % default for bounded BCI is not to keep params across trials
end

digitalCodeNameBciStartsAfter = expParams.bciStartsAfterCode;
digitalCodeTrialStart = codes.(digitalCodeNameBciStartsAfter);% could be START_TRIAL
digitalCodeNameBciEndsBy = expParams.bciEndsByCode;
digitalCodeTrialEnd = codes.(digitalCodeNameBciEndsBy);% could be END_TRIAL

% this could be the code associated with FIX_OFF or TARG_OFF for example,
% or it could be a BCI specific code. Leaving it open which would be
% best... these are distinct from trial start/end, because *those* allow us
% to explore stuff in the intertrial period, while these tell us to focus
% on the BCI
digitalCodeNameBciStart = expParams.bciStartCode;
digitalCodeBciStart = codes.(digitalCodeNameBciStart);
digitalCodeNameBciEnd = expParams.bciEndCode;
digitalCodeBciEnd = codes.(digitalCodeNameBciEnd);

% Determine whether online mat containing received spikes should be saved
% or not.
if saveOnlineMatFile
    % Keep track of BCI_Correct and BCI_Missed trials; needed for offline analysis 
    digitalCodeNameBciCorrect = expParams.bciCorrectCode;
    digitalCodeNameBciMissed = expParams.bciMissedCode;
    digitalCodeNameBciAbort= expParams.bciAbortCode;

    %digitalCodeNameBciIncorrect = expParams.bciIncorrectCode;
    digitalCodeBciCorrect = codes.(digitalCodeNameBciCorrect);
    digitalCodeBciMissed = codes.(digitalCodeNameBciMissed);
    digitalCodeBciAbort = codes.(digitalCodeNameBciAbort);
    
    taskName = expParams.exFileName;
    onlineBCIMatStruct = struct('trialIdx', {}, 'trialSpikes', {}, 'trialReturnVals', {}, 'trialType', {}, 'bciTrialResult', {}, 'bciTrialParams', {});
    trialIdx = 0;
    onlineMatFileName= sprintf('%s_%s_%sOnlineDat.mat', datestr(today, 'yyyymmdd'), datestr(now, 'HH-MM-SS'), taskName);
    onlineMatDir = sprintf('%s/%s',params.bciDecoderBasePathBciComputer, expParams.subject);
    if not(isfolder(onlineMatDir))
        mkdir(onlineMatDir)
    end
    fullFileName = fullfile(onlineMatDir, onlineMatFileName);
    fprintf('Online File saved at : %s \n', fullFileName)
end

boundStarted = false;
customBciCodeAfterTrlStart = [];
bciTrialResult = [];
timePtBoundStarted = [];
timePtBoundEnded = [];
bciStart = false;
timePtBciStarted = [];
timePtBciEnd = [];
timePtBciResult = [];

samplesPerSecond = params.neuralRecordingSamplingFrequencyHz;%30000;
binSizeMs = expParams.binSizeMs;%50;
nasNetwork = expParams.nasNetwork;
currReturn = expParams.initReturn'; % i.e. [0,0] if velocity
[nasNetParams.w1, nasNetParams.b1, nasNetParams.w2, nasNetParams.b2] = loadNasNet(nasNetwork);
gamma = expParams.gamma;

bciDecoderFunctionName = expParams.name;
bciDecoderFunction = str2func(bciDecoderFunctionName);
if refreshOutput, clear(bciDecoderFunctionName); end % make sure it's fresh

msPerS = 1000;
samplesPerBin = binSizeMs/msPerS*samplesPerSecond;

binSpikeCountOverall = zeros(length(okelecs),1);
binSpikeCountNextOverall = zeros(length(okelecs), 1);
delValues = '';

% DEBUGGING
binCntNasTrial = {};
allTmstmpAll = {};
tmstpInit = [];
waveforms = [];
binNum = -1;

% grab events
% [count,tmstp,events]=xippmex('digin');
% prlEvents = [events.parallel];
modelParams = [];
loopTmTotalSec = 0;
while true
    loopTmStart = tic;
    % check for messages or BCI end between trials...
    [bciEnd, ctrlMsg] = checkIfBciEndOrMsg(controlCompSocket);
    % check the start of thectrlMsg trial using codes sent to Ripple
    [count,tmstpPrlEvt,events]=xippmex('digin');
    % Is a buffer where we get some codes at a time (Not continuous access
    % to whole Nev)
    prlEvents = [events.parallel];
    tstpTrlStart = find(prlEvents==digitalCodeTrialStart);
    tstpTrlEnd = find(prlEvents==digitalCodeTrialEnd);
    tstpBciStart = find(prlEvents==digitalCodeBciStart);
    tstpBciEnd = find(prlEvents==digitalCodeBciEnd);
    if saveOnlineMatFile
        % Find indices of prlEvents at which bciCorrect or bciIncorrect code is
        % sent.
        tstpBciCorrect = find(prlEvents == digitalCodeBciCorrect);
        tstpBciMissed = find(prlEvents == digitalCodeBciMissed);
        tstpBciAbort = find(prlEvents == digitalCodeBciAbort);
    end
    % Find indices at which custom codes were sent
    tstpCustomBciCode = find((prlEvents>20000) & (prlEvents<20010));
    % Only keep custom BCI Codes that occur after current trial's start
    % timestamp
    customBciCode = prlEvents(((prlEvents>20000) & (prlEvents<20010)));
    if length(tstpTrlStart)>1
        disp('missed a trial')
        tstpTrlStart = tstpTrlStart(end);
        % Force increment trial counter if missed 
        trialIdx = trialIdx + 1; 
    end
    if ~isempty(tstpTrlStart)
        disp('trial start')
        timePtBoundStarted = tmstpPrlEvt(tstpTrlStart);
        boundStarted = true;
         % Initialize arrays that we'll  track during trials
        currTrialSpikesArray = [];
        currTrialReturnVals = [];
        % Initialize trial type and result to be 0
        currTrialTypeIdx = 0;
        currTrialResult = 0;
        currTrialBCIParams = [];
    end
    
    % Set the custom BCI Code
    if(~isempty(tstpCustomBciCode))
        % Find actual time points at which customBCICodes were sent
        timePtCustomBciCode = tmstpPrlEvt(tstpCustomBciCode);
        % Only pick up custom indices After the bound has started
        indicesForCustomCodesSentAfterBound = timePtCustomBciCode>timePtBoundStarted;
        % Find Timestamps of custom bci code that are after bound start
        tstpCustomBciCodeAfterBoundStart = tstpCustomBciCode(indicesForCustomCodesSentAfterBound);
        if ~isempty(tstpCustomBciCodeAfterBoundStart)
            % Set it to the first index if multiple codes sent for some
            % reason
            customBCICodeIdx = tstpCustomBciCodeAfterBoundStart(1);
            % Make sure that customBCICodeIdx doesn't contain
            % indices outside of customBCICode
            if ~isempty(customBciCode) && (length(customBciCode) >= customBCICodeIdx)
                customBciCodeAfterTrlStart = customBciCode(customBCICodeIdx)-20000;
            end
        end
    end
    if saveOnlineMatFile
        % Added by Chris 11-20-23
    %     % Check if either BCI correct or bci Missed is received
        if(~isempty(tstpBciCorrect))
            % Find actual time points at which BCI result was sent
            timePtBciResult = tmstpPrlEvt(tstpBciCorrect);
            % Only pick up result indices After the bound has started
            indicesForResultsSentAfterBound = timePtBciResult>timePtBoundStarted;
            if(sum(indicesForResultsSentAfterBound)>0)
                currTrialResult=digitalCodeBciCorrect;
            end
        elseif(~isempty(tstpBciMissed))
            timePtBciResult = tmstpPrlEvt(tstpBciMissed);
            % Only pick up custom indices After the bound has started
            indicesForResultsSentAfterBound = timePtBciResult>timePtBoundStarted;
            if(sum(indicesForResultsSentAfterBound)>0)
                currTrialResult=digitalCodeBciMissed;
            end
        elseif(~isempty(tstpBciAbort))
            timePtBciResult = tmstpPrlEvt(tstpBciAbort);
            % Only pick up custom indices After the bound has started
            indicesForResultsSentAfterBound = timePtBciResult>timePtBoundStarted;
            if(sum(indicesForResultsSentAfterBound)>0)
                currTrialResult=digitalCodeBciAbort;
            end
        end
    end
    
    if length(tstpTrlEnd)>1
        % might happen with length(tstpTrlStart)==1 if we catch the end of
        % the previous and the start and end of the current
        tstpTrlStart = tstpTrlStart(end);
    end
    if ~isempty(tstpTrlEnd)
        timePtBoundEnded = tmstpPrlEvt(tstpTrlEnd);
    end
    
    if timePtBoundEnded > timePtBoundStarted
        if boundStarted
            disp('trial end')
            % Save copy of online mat at end of every trial
            if saveOnlineMatFile
                % Append to online BCI struct array the current trial's
                % information
                onlineBCIMatStruct(end+1) = struct('trialIdx', trialIdx, 'trialSpikes', currTrialSpikesArray , 'trialReturnVals', currTrialReturnVals, 'trialType', currTrialTypeIdx, 'bciTrialResult', currTrialResult,'bciTrialParams', currTrialBCIParams);
                save(fullFileName, 'onlineBCIMatStruct', '-v6');
                % Increment trial counter at end
                trialIdx = trialIdx + 1; 
            end
            
        end
        if bciStart
            % DEBUGGING
%             binCntNasTrial = [binCntNasTrial binSpkCntTrial];
%             allTmstmpAll = [allTmstmpAll {allTmstmpTrl}];
            % END DEBUGGING

            fprintf('bci end  with trial after %d bins\n', binNum)
            binNum = -1;
        end
        boundStarted = false;
        customBciCodeAfterTrlStart = [];
        bciStart = false;
        currReturn = expParams.initReturn';
        if refreshOutput, clear(bciDecoderFunctionName); end % in a bounded BCI, we clear persistent variables after the end of the bound
    end
  
    if boundStarted      
        [modelParams, updatedReturn] = processBciControlMessage(controlCompSocket, ctrlMsg, modelParams);
        if ~isempty(updatedReturn)
            currReturn = updatedReturn;
        end
        
        % buffering issues cause weird timing--specifically, some channels
        % will have smaller timestamps then the previous call of other
        % channels; I think the smaller the buffer the less this is a
        % problem, so for now calling this every cycle I think is a good
        % idea
        prevTmstpInit = tmstpInit;
        prevWaveforms = waveforms;
        [~,tmstpInit, waveforms, ~]=xippmex('spike',okelecs,zeros(1,length(okelecs)));
    
        if ~isempty(modelParams)
            % in case we have two starts/ends, we only want the start related to the current trial
                        
            if ~isempty(tstpBciStart)
                disp('bci start')
                timePtBciStarted = tmstpPrlEvt(tstpBciStart);
                timePtBciStarted = timePtBciStarted(timePtBciStarted>=timePtBoundStarted);
                timePtBinStart = timePtBciStarted;
                goodChannelNums = modelParams.channelsKeep;
                goodChannelInds = ismember(okelecs, goodChannelNums);
                binSpikeCountOverall = zeros(length(goodChannelNums), 1);
                binSpikeCountNextOverall = zeros(length(goodChannelNums), 1);
                bciStart = true;
                bciJustStarted = true; % important for grabbing any early spikes
                % Initialize arrays that we'll  track during BCI trials
                currTrialSpikesArray = zeros(length(goodChannelNums), 0);
                currTrialReturnVals = zeros(length(currReturn), 0);
                currTrialResult = 0;
                % Trial Type set to zero if not initialized
                currTrialTypeIdx = 0;
                % Set BCITrialParams to -99999
                currTrialBCIParams = [];
            end
            if ~isempty(tstpBciEnd)
                timePtBciEnd = tmstpPrlEvt(tstpBciEnd);
                timePtBciEnd = timePtBciEnd(timePtBciEnd>timePtBoundStarted);
            end
            
            if bciStart  
                if bciJustStarted
                    binNum=0;
                    [countsPerChannel, countsPerChannelNextBin] = countBinnedSpikesPerChannel(prevTmstpInit, goodChannelInds, timePtBinStart, samplesPerBin, nasNetParams, prevWaveforms, gamma);
                    binSpikeCountOverall = binSpikeCountOverall + countsPerChannel;
                    % this'll likely just be zero most of the time, but on bin
                    % edges it'll be needed
                    binSpikeCountNextOverall = binSpikeCountNextOverall + countsPerChannelNextBin;
                    bciJustStarted = false;
                end
                
                % DEBUGGING
%                 allTmstmpTrl = [allTmstmpTrl tmstpInit(goodChannelInds)];
                % END DEBUGGING
                [countsPerChannel, countsPerChannelNextBin] = countBinnedSpikesPerChannel(tmstpInit, goodChannelInds, timePtBinStart, samplesPerBin, nasNetParams, waveforms, gamma);
                binSpikeCountOverall = binSpikeCountOverall + countsPerChannel;
                % this'll likely just be zero most of the time, but on bin
                % edges it'll be needed
                binSpikeCountNextOverall = binSpikeCountNextOverall + countsPerChannelNextBin;
                
                % allTmstps is being used to check bin turnover, rather
                % than counting spikes, so we want to see as many
                % timestamps as possible here, so as not to miss a bin
                % (bonus that BCI_CURSOR_POS shoutimePtBciEndld always be getting sent
                % within 50ms, so tmstpPrlEvt is definitely going to give
                % good bin cutoffs)
                allTmstps = cat(2, tmstpInit{:}, tmstpPrlEvt);
                if any(allTmstps>(timePtBinStart+2*samplesPerBin))
                    % this might happen if the recorded waveforms include
                    % the next bin, and also have waveforms from two bins
                    % after, but we're noting we should re    ally expect at
                    % most samples in the next bin
                    fprintf('furthest out sample (shooould be less than %d at most) is %d\n', 2*samplesPerBin, max(allTmstps)-timePtBinStart)
                end
                
                % by checking whether there are timestamps from the *next*
                % bin, we confirm that we've completed the current bin and
                % can send off info
                if any(allTmstps>(timePtBinStart+samplesPerBin))
                    if ~any(allTmstps<=(timePtBinStart+samplesPerBin))
                        fprintf('max (samples, time) past bin #%d end: (%d, %0.2f ms)\n', binNum, max(allTmstps-(timePtBinStart+samplesPerBin)), max(allTmstps-(timePtBinStart+samplesPerBin))/samplesPerSecond*msPerS);
                    end
                    
                    % but actually xippmex kind of lies because of how
                    % threshold crossings are built up--it cycles through
                    % each channel and updates the buffer, so sometimes one
                    % channel has spikes from the next bin but another
                    % channel hasn't quite buffered its own spikes from the
                    % *current* bin. Hopefully by calling xippmex one more
                    % time we can catch those extra spikes and make it
                    % vanishingly unlikely that we miss something...
                    [~,tmstpInit, waveforms, ~]=xippmex('spike',okelecs,zeros(1,length(okelecs)));
                    % DEBUGGING
%                     allTmstmpTrl = [allTmstmpTrl tmstpIn    it(goodChannelInds)];
                    % END DEBUGGING
                    [countsPerChannel, countsPerChannelNextBin] = countBinnedSpikesPerChannel(tmstpInit, goodChannelInds, timePtBinStart, samplesPerBin, nasNetParams, waveforms, gamma);
                    % scoop up last current bin spikes and grow the next
                    % bin spikes
                    binSpikeCountOverall = binSpikeCountOverall + countsPerChannel;
                    binSpikeCountNextOverall = binSpikeCountNextOverall + countsPerChannelNextBin;
                    
                    meanSpikeCount = mean(binSpikeCountOverall,2);
                    
                    % run the BCI decoder
                    currReturnSend = [currReturn; customBciCodeAfterTrlStart];
                    currReturn = bciDecoderFunction(meanSpikeCount, currReturnSend, modelParams, expParams, binNum, controlCompSocket);
                    
                    % Keep track of current spike counts and currReturn
                    if saveOnlineMatFile
                        currTrialSpikesArray(:, end+1) = meanSpikeCount;
                        currTrialReturnVals(:,end+1) = currReturn;
                        % Specific to reward axis BCI (make this more
                        % generalizable)
                        currTrialBCIParams = [expParams.smallTargChangeByStd, expParams.largeTargChangeByStd];
                        currTrialTypeIdx = customBciCodeAfterTrlStart;
                    end
                    % prep the message to send
                    uint8Msg = typecast(currReturn, 'uint8');
                    if size(uint8Msg, 1) ~= 1
                        msgToSend = uint8Msg';
                    else
                        msgToSend = uint8Msg;
                    end
                    matlabUDP2('send',controlCompSocket.sender, msgToSend);

                    % Keep track of currMeanSpikes and currReturn
                    
                    
                    % the current bin is now what was the next bin before
                    binSpikeCountOverall = binSpikeCountNextOverall;

                    % zero out the counts for the next bin
                    binSpikeCountNextOverall(:) = 0;
                    timePtBinStart = timePtBinStart+samplesPerBin;
                    binNum = binNum+1; % keep track of bin number
                end
            end
            
            % allow the BCI loop to run one final time to see if BCI ended
            % after a full bin happened
            if timePtBciEnd > timePtBciStarted
                if bciStart
                    fprintf('bci end in trial after %d bins\n', binNum)
                    binNum = -1;
                end
                bciStart = false;
                currReturn = expParams.initReturn';
                if refreshOutput, clear(bciDecoderFunctionName); end % in a bounded BCI, we clear persistent variables after the end of the bound
            end
        end
    end
    loopTmTotalSec = toc(loopTmStart);
    if loopTmTotalSec>binSizeMs/1000
        fprintf('loop time of %d in bin number %d longer than binSizeMs of %d\n', loopTmTotalSec, binNum, binSizeMs);
        if binNum<0
            fprintf('bin number of -1 means this happened outside the BCI\n');
        end
    end
end
end
