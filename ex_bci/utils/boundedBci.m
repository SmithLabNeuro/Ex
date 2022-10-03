function boundedBci(controlCompSocket, expParams, okelecs)

global params codes
global tmstpNasSpkAll
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

boundStarted = false;
timePtBoundStarted = [];
timePtBoundEnded = [];
bciStart = false;
timePtBciStarted = [];
timePtBciEnd = [];
samplesPerSecond = params.neuralRecordingSamplingFrequencyHz;%30000;
binSizeMs = expParams.binSizeMs;%50;
nasNetwork = expParams.nasNetwork;
currReturn = expParams.initReturn'; % i.e. [0,0] if velocity
[nasNetParams.w1, nasNetParams.b1, nasNetParams.w2, nasNetParams.b2] = loadNasNet(nasNetwork);
gamma = expParams.gamma;

bciDecoderFunctionName = expParams.name;
bciDecoderFunction = str2func(bciDecoderFunctionName);
clear(bciDecoderFunctionName); % make sure it's fresh

msPerS = 1000;
samplesPerBin = binSizeMs/msPerS*samplesPerSecond;

binSpikeCountOverall = zeros(length(okelecs),1);
binSpikeCountNextOverall = zeros(length(okelecs), 1);
delValues = '';

tmstpNasSpkAll = [];
binCntNasTrial = {};
allTmstmpAll = {};
tmstpInit = [];
waveforms = [];
binNum = 0;

% grab events
% [count,tmstp,events]=xippmex('digin');
% prlEvents = [events.parallel];
modelParams = [];
loopTmAll = 0;
while true
    
    loopTm = tic;
    % check for messages or BCI end between trials...
    [bciEnd, ctrlMsg] = checkIfBciEndOrMsg(controlCompSocket);
    if bciEnd
        break
    end
    % check the start of the trial using codes sent to Ripple
    [count,tmstpPrlEvt,events]=xippmex('digin');
    prlEvents = [events.parallel];
    tstpTrlStart = find(prlEvents==digitalCodeTrialStart);
    tstpTrlEnd = find(prlEvents==digitalCodeTrialEnd);
    tstpBciStart = find(prlEvents==digitalCodeBciStart);
    tstpBciEnd = find(prlEvents==digitalCodeBciEnd);
    
    if length(tstpTrlStart)>1
        disp('missed a trial')
        tstpTrlStart = tstpTrlStart(end);
    end
    if ~isempty(tstpTrlStart)
        disp('trial start')
        timePtBoundStarted = tmstpPrlEvt(tstpTrlStart);
        boundStarted = true;
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
%             save('/home/smithlab/tempChecker.mat', 'binCntNasTrial', 'allTmstmpAll');
            disp('trial end')
        end
        if bciStart
            binCntNasTrial = [binCntNasTrial binSpkCntTrial];
            allTmstmpAll = [allTmstmpAll {allTmstmpTrl}];

            fprintf('bci end  with trial after %d bins\n', binNum)
        end
        boundStarted = false;
        bciStart = false;
        currReturn = expParams.initReturn';
        clear(bciDecoderFunctionName); % in a bounded BCI, we clear persistent variables after the end of the bound
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
                % DEBUGGING
                binSpkCntTrial = zeros(length(goodChannelNums), 0);
                allTmstmpTrl = [];
            end
            if ~isempty(tstpBciEnd)
                timePtBciEnd = tmstpPrlEvt(tstpBciEnd);
                timePtBciEnd = timePtBciEnd(timePtBciEnd>timePtBoundStarted);
            end
            
            if bciStart  
                if bciJustStarted
                    binNum=1;
                    [countsPerChannel, countsPerChannelNextBin] = countBinnedSpikesPerChannel(prevTmstpInit, goodChannelInds, timePtBinStart, samplesPerBin, nasNetParams, prevWaveforms, gamma);
                    binSpikeCountOverall = binSpikeCountOverall + countsPerChannel;
                    % this'll likely just be zero most of the time, but on bin
                    % edges it'll be needed
                    binSpikeCountNextOverall = binSpikeCountNextOverall + countsPerChannelNextBin;
                    bciJustStarted = false;
                else
                    binNum = binNum+1;
                end
                
                allTmstmpTrl = [allTmstmpTrl tmstpInit(goodChannelInds)];
                [countsPerChannel, countsPerChannelNextBin] = countBinnedSpikesPerChannel(tmstpInit, goodChannelInds, timePtBinStart, samplesPerBin, nasNetParams, waveforms, gamma);
                binSpikeCountOverall = binSpikeCountOverall + countsPerChannel;
                % this'll likely just be zero most of the time, but on bin
                % edges it'll be needed
                binSpikeCountNextOverall = binSpikeCountNextOverall + countsPerChannelNextBin;
                
                % allTmstps is being used to check bin turnover, rather
                % than counting spikes, so we want to see as many
                % timestamps as possible here, so as not to miss a bin
                % (bonus that BCI_CURSOR_POS should always be getting sent
                % within 50ms, so tmstpPrlEvt is definitely going to give
                % good bin cutoffs)
                allTmstps = cat(2, tmstpInit{:}, tmstpPrlEvt);
                if any(allTmstps>(timePtBinStart+2*samplesPerBin))
                    % this might happen if the recorded waveforms include
                    % the next bin, and also have waveforms from two bins
                    % after, but we're noting we should really expect at
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
                    allTmstmpTrl = [allTmstmpTrl tmstpInit(goodChannelInds)];
                    [countsPerChannel, countsPerChannelNextBin] = countBinnedSpikesPerChannel(tmstpInit, goodChannelInds, timePtBinStart, samplesPerBin, nasNetParams, waveforms, gamma);
                    % scoop up last current bin spikes and grow the next
                    % bin spikes
                    binSpikeCountOverall = binSpikeCountOverall + countsPerChannel;
                    binSpikeCountNextOverall = binSpikeCountNextOverall + countsPerChannelNextBin;
                    
                    meanSpikeCount = mean(binSpikeCountOverall,2);
                    % DEBUGGING
                    binSpkCntTrial = [binSpkCntTrial meanSpikeCount];
                    tmstpNasSpkAll = [];
                    % END DEBUGGING
                    
                    % run the BCI decoder
                    currReturn = bciDecoderFunction(meanSpikeCount, currReturn, modelParams, expParams);
                    
                    % prep the message to send
                    uint8Msg = typecast(currReturn, 'uint8');
                    if size(uint8Msg, 1) ~= 1
                        msgToSend = uint8Msg';
                    else
                        msgToSend = uint8Msg;
                    end
                    matlabUDP2('send',controlCompSocket.sender, msgToSend);

                    % the current bin is now what was the next bin before
                    binSpikeCountOverall = binSpikeCountNextOverall;
                    
                    % zero out the counts for the next bin
                    binSpikeCountNextOverall(:) = 0;
                    timePtBinStart = timePtBinStart+samplesPerBin;
                end
            end
            
            % allow the BCI loop to run one final time to see if BCI ended
            % after a full bin happened
            if timePtBciEnd > timePtBciStarted
                if bciStart
                    binCntNasTrial = [binCntNasTrial binSpkCntTrial];
                    allTmstmpAll = [allTmstmpAll {allTmstmpTrl}];
                    fprintf('bci end in trial after %d bins\n', binNum)
                end
                bciStart = false;
                currReturn = expParams.initReturn';
                clear(bciDecoderFunctionName); % in a bounded BCI, we clear persistent variables after the end of the bound
            end
        end
    end
    loopTmAll = toc(loopTm);
end

save('/home/smithlab/tempChecker.mat', 'binCntNasTrial', 'allTmstmpAll');