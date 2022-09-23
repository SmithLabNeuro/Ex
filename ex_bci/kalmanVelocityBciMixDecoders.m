function kalmanVelocityDualDecoderBci(controlCompSocket, expParams, okelecs)

global params codes

digitalCodeTrialStart = codes.START_TRIAL;%1;
digitalCodeTrialEnd = codes.END_TRIAL;%255;

% this could be the code associated with FIX_OFF or TARG_OFF for example,
% or it could be a BCI specific code. Leaving it open which would be
% best... these are distinct from trial start/end, because *those* allow us
% to explore stuff in the intertrial period, while these tell us to focus
% on the BCI
digitalCodeNameBciStart = expParams.bciStartCode;
digitalCodeBciStart = codes.(digitalCodeNameBciStart);
digitalCodeNameBciEnd = expParams.bciEndCode;
digitalCodeBciEnd = codes.(digitalCodeNameBciEnd);

trialStarted = false; 
timePtStarted = [];
timePtEnded = [];
bciStart = false;
timePtBciStarted = [];
timePtBciEnd = [];
samplesPerSecond = params.neuralRecordingSamplingFrequencyHz;%30000;
binSizeMs = expParams.binSizeMs;%50;
nasNetwork = expParams.nasNetwork;
weightedVelocity = [0; 0];
[w1, b1, w2, b2] = loadNasNet(nasNetwork);
firstVelScale = expParams.firstVelScale;
secondVelScale = expParams.secondVelScale;
gamma_1 = expParams.firstGamma;
gamma_2 = expParams.secondGamma;
decoderParameterLocation = params.bciDecoderBasePathBciComputer;

msPerS = 1000;
samplesPerBin = binSizeMs/msPerS*samplesPerSecond;

binSpikeCountOverall_1 = zeros(length(okelecs),1);
binSpikeCountOverall_2 = zeros(length(okelecs),1);
delValues = '';

% grab events
% [count,tmstp,events]=xippmex('digin');
% prlEvents = [events.parallel];
modelParams_1 = [];
modelParams_2 = [];
loopTmAll = 0;
scaleIndex = 1;
while true
    
    loopTm = tic;
    % check for messages or BCI end between trials...
    [bciEnd, ctrlMsg] = checkIfBciEndOrMsg(controlCompSocket);
    if bciEnd
        break
    end
    % check the start of the trial using codes sent to Ripple
    a = tic;
    [count,tmstpPrlEvt,events]=xippmex('digin');
    prlEvents = [events.parallel];
    tstpTrlStart = find(prlEvents==digitalCodeTrialStart);
    tstpTrlEnd = find(prlEvents==digitalCodeTrialEnd);
    tstpBciStart = find(prlEvents==digitalCodeBciStart);
    tstpBciEnd = find(prlEvents==digitalCodeBciEnd);
    prlTm = toc(a);
    
    if length(tstpTrlStart)>1
        disp('missed a trial')
        tstpTrlStart = tstpTrlStart(end);
    end 
    if ~isempty(tstpTrlStart)
        disp('trial start')
        timePtStarted = tmstpPrlEvt(tstpTrlStart);
        trialStarted = true;
    end
    
    if length(tstpTrlEnd)>1
        % might happen with length(tstpTrlStart)==1 if we catch the end of
        % the previous and the start and end of the current
        tstpTrlStart = tstpTrlStart(end);
    end
    if ~isempty(tstpTrlEnd)
        timePtEnded = tmstpPrlEvt(tstpTrlEnd);
    end
    
    if timePtEnded > timePtStarted
        if trialStarted
            disp('trial end')
        end
        if bciStart
            disp('bci end with trial')
        end
        trialStarted = false;
        bciStart = false;
        weightedVelocity = [0; 0];
    end
  
    if trialStarted
        if ~isempty(ctrlMsg)
            % check if there's a new decoder or if there's an adjusted
            % velocity! (Important for calibration steps!)
            if strcmp(ctrlMsg, 'decoderParameterFile')
                decoderParameterFileRelativePath_1 = receiveMessageSendAck(controlCompSocket);
                decoderParameterFileRelativePath_1(decoderParameterFileRelativePath_1=='\') = '/';
%                 decoderParameterFileRelativePath = 'satchel/21-Mar-2022/KalmanBci_14-02-09.mat';
                decoderParameterFileFullPath_1 = fullfile(decoderParameterLocation, decoderParameterFileRelativePath_1);
                modelParams_1 = load(decoderParameterFileFullPath_1, 'M0', 'M1', 'M2', 'channelsKeep');
                M0_1 = modelParams_1.M0;
                M1_1 = modelParams_1.M1;
                M2_1 = modelParams_1.M2;
                goodChannelNums_1 = modelParams_1.channelsKeep;
                goodChannelInds_1 = ismember(okelecs, goodChannelNums_1);
                fprintf('loaded new parameters from %s\n', decoderParameterFileRelativePath_1)

                decoderParameterFileRelativePath_2 = receiveMessageSendAck(controlSocket);
                decoderParameterFileRelativePath_2(decoderParameterFileRelativePath_2=='\') = '/';
                decoderParameterFileFullPath_2 = fullfile(decoderParameterLocation, decoderParameterFileRelativePath_2);
                modelParams_2 = load(decoderParameterFileFullPath_2, 'M0', 'M1', 'M2', 'channelsKeep');                
                M0_2 = modelParams_2.M0;
                M1_2 = modelParams_2.M1;
                M2_2 = modelParams_2.M2; 
                goodChannelNums_2 = modelParams_2.channelsKeep;
                goodChannelInds_2 = ismember(okelecs, goodChannelNums_2);
                fprintf('loaded new parameters from %s\n', decoderParameterFileRelativePath_2)
            elseif strcmp(ctrlMsg, 'scaleIndex')
                scaleIndex = str2double(receiveMessageSendAck(controlCompSocket));  
                fprintf('NEW SCALE INDEX: %d\n', scaleIndex);
            else
                weightedVelocity = typecast(uint8(ctrlMsg), 'double')';
%                 fprintf('received constrained velocity\n');
            end
        end
        if ~isempty(modelParams_1) && ~isempty(modelParams_2)
            % in case we have two starts/ends, we only want the start related to the current trial
                        
            if ~isempty(tstpBciStart)
                disp('bci start')
                timePtBciStarted = tmstpPrlEvt(tstpBciStart);
                timePtBciStarted = timePtBciStarted(timePtBciStarted>timePtStarted);
                timePtBinStart = timePtBciStarted;
                binSpikeCountOverall_1 = zeros(length(goodChannelNums_1), 1);
                binSpikeCountOverall_2 = zeros(length(goodChannelNums_2), 1);
                bciStart = true;
            end
            if ~isempty(tstpBciEnd)
                timePtBciEnd = tmstpPrlEvt(tstpBciEnd);
                timePtBciEnd = timePtBciEnd(timePtBciEnd>timePtStarted);
            end
            % Checking to make sure that end of BCI trial is always
            % received before the BCI start trial code.
            if timePtBciEnd > timePtBciStarted
                if bciStart
                    disp('bci end in trial')
                end
                bciStart = false;
                weightedVelocity = [0; 0];
            end
            
            if bciStart
                % Note to self: would it work to set okelecs to
                % goodChannels in some way?
                [~,tmstpInit, waveforms, ~]=xippmex('spike',okelecs,zeros(1,length(okelecs)));
                
                
                % only use data from channels that were good for
                % calibration
                tmstpGoodChSpk_1 = tmstpInit(goodChannelInds_1);
                waveforms_1 = waveforms(goodChannelInds_1);
                
                tmstpGoodChSpk_2 = tmstpInit(goodChannelInds_2);
                waveforms_2 = waveforms(goodChannelInds_2);
                
                % run waveforms through NAS net
                tmstpNasSpk_1 = cellfun(@(wvForms, tms) tms(runNASNetContinuous(w1, b1, w2, b2, wvForms, gamma_1)), waveforms_1, tmstpGoodChSpk_1, 'uni', 0);
%                 tmstp = cellfun(@(wvForms, tms) tms, waveforms, tmstp, 'uni', 0);
                spikesThisBinByChannel_1 = cellfun(@(x) x>timePtBinStart & x<timePtBinStart+samplesPerBin, tmstpNasSpk_1, 'uni', 0);
                %                 waveformsThisBinByChannel = cellfun(@(wvFrm, spksInBin) wvFrm(spksInBin, :), waveforms, spikesThisBinByChannel, 'uni', 0);
                spikesNextBinByChannel_1 = cellfun(@(x) x>timePtBinStart+samplesPerBin, tmstpNasSpk_1, 'uni', 0);
                %                 waveformsNextBinByChannel = cellfun(@(wvFrm, spksInBin) wvFrm(spksInBin, :), waveforms, spikesNextBinByChannel, 'uni', 0);
                countsPerChannelCell_1 = cellfun(@(x) sum(x, 2), spikesThisBinByChannel_1, 'uni', 0);
                countsExistChannel_1 = ~cellfun('isempty', countsPerChannelCell_1);
                countsPerChannel_1 = zeros(length(countsPerChannelCell_1), 1);
                countsPerChannel_1(countsExistChannel_1) = [countsPerChannelCell_1{countsExistChannel_1}];
                binSpikeCountOverall_1 = binSpikeCountOverall_1 + countsPerChannel_1;

                
                % Do same thing but for second set of waveforms
                tmstpNasSpk_2 = cellfun(@(wvForms, tms) tms(runNASNetContinuous(w1, b1, w2, b2, wvForms, gamma_2)), waveforms_2, tmstpGoodChSpk_2, 'uni', 0);
%                 tmstp = cellfun(@(wvForms, tms) tms, waveforms, tmstp, 'uni', 0);
                spikesThisBinByChannel_2 = cellfun(@(x) x>timePtBinStart & x<timePtBinStart+samplesPerBin, tmstpNasSpk_2, 'uni', 0);
                %                 waveformsThisBinByChannel = cellfun(@(wvFrm, spksInBin) wvFrm(spksInBin, :), waveforms, spikesThisBinByChannel, 'uni', 0);
                spikesNextBinByChannel_2 = cellfun(@(x) x>timePtBinStart+samplesPerBin, tmstpNasSpk_2, 'uni', 0);
                %                 waveformsNextBinByChannel = cellfun(@(wvFrm, spksInBin) wvFrm(spksInBin, :), waveforms, spikesNextBinByChannel, 'uni', 0);
                countsPerChannelCell_2 = cellfun(@(x) sum(x, 2), spikesThisBinByChannel_2, 'uni', 0);
                countsExistChannel_2 = ~cellfun('isempty', countsPerChannelCell_2);
                countsPerChannel_2 = zeros(length(countsPerChannelCell_2), 1);
                countsPerChannel_2(countsExistChannel_2) = [countsPerChannelCell_2{countsExistChannel_2}];
                binSpikeCountOverall_2 = binSpikeCountOverall_2 + countsPerChannel_2;
                
                % allTmstps is being used to check bin turnover, rather
                % than counting spikes, so we want to see as many
                % timestamps as possible here, so as not to miss a bin
                % (bonus that BCI_CURSOR_POS should always be getting sent
                % within 50ms, so tmstpPrlEvt is definitely going to give
                % good bin cutoffs) %
                % tmstpInit represents all waveforms (not channel filtered) filtered by NAS and
                % tmstpPrlEvent represents all digital codes 
                allTmstps = cat(2, tmstpInit{:}, tmstpPrlEvt);
                if any(allTmstps>(timePtBinStart+2*samplesPerBin))
                    if any(allTmstps<=(timePtBinStart+samplesPerBin))
                        % this might happen if the recorded waveforms come from this
                        % bin, include the next bin, and also have waveforms from two
                        % bins after...
                        fprintf('furthest out sample (shooould be less than %d) is %d\n', samplesPerBin, max(allTmstps)-timePtBinStart)
                        %                     warning('BCI code is running too slow for bin time...')
                        
                        %                     fprintf('prl read ran in %0.3f secs\n', prlTm);
                        %                     fprintf('loop ran in %0.3f secs\n', loopTmAll);
%                     fprintf('loop current is taking %0.3f secs\n', toc(loopTm));

                    else
                        fprintf('why am I here?\n');
%                         timePtBinStart = timePtBinStart + samplesPerBin * floor(max(allTmstps - timePtBinStart)/samplesPerBin);
                    end
                end
                if any(allTmstps>(timePtBinStart+samplesPerBin))
                    if ~any(allTmstps<=(timePtBinStart+samplesPerBin))
                        fprintf('samples past bin end: %d\n', min(allTmstps-(timePtBinStart+samplesPerBin)));
                    end
                    % by checking whether there are timestamps from the
                    % *next* bin, we confirm that we've completed the
                    % current bin and can send off info
                    meanSpikeCount_1 = mean(binSpikeCountOverall_1,2);
                    meanSpikeCount_2 = mean(binSpikeCountOverall_2,2);
                    
                    % Compute velocities separately for both decoders
                    velocity_1 = M0_1 + M1_1*weightedVelocity + M2_1 * meanSpikeCount_1;
                    velocity_2 = M0_2 + M1_2*weightedVelocity + M2_2* meanSpikeCount_2;
                    weightedVelocity = firstVelScale(scaleIndex)*velocity_1 + secondVelScale(scaleIndex)*velocity_2;
                    uint8Msg = typecast(weightedVelocity, 'uint8');
                    if any(uint8Msg==0)
                        %                     disp(uint8Msg)
                        %                     disp(typecast(meanSpikeCount, 'int8'))
                    end
                    msgToSend = uint8Msg';
%                     disp(velocity')
                    
                    matlabUDP2('send',controlCompSocket.sender, msgToSend);
%                     fprintf('sent unconstrained velocity\n');
                    
                    % reset the bin spikes by computing how much of the next bin
                    % was captured here, and shift the starting point
                    countsPerChannelNextBin_1 = zeros(size(binSpikeCountOverall_1));
                    emptyNextBin_1 = cellfun('isempty', spikesNextBinByChannel_1);
                    countsPerChannelNextBin_1(~emptyNextBin_1) = cellfun(@(x) sum(x, 2), spikesNextBinByChannel_1(~emptyNextBin_1));
                    binSpikeCountOverall_1 = countsPerChannelNextBin_1;
                    timePtBinStart = timePtBinStart+samplesPerBin;
                    
                    countsPerChannelNextBin_2 = zeros(size(binSpikeCountOverall_2));
                    emptyNextBin_2 = cellfun('isempty', spikesNextBinByChannel_2);
                    countsPerChannelNextBin_2(~emptyNextBin_2) = cellfun(@(x) sum(x, 2), spikesNextBinByChannel_2(~emptyNextBin_2));
                    binSpikeCountOverall_2 = countsPerChannelNextBin_2;
                end
            end
        end
    end
    loopTmAll = toc(loopTm);
end