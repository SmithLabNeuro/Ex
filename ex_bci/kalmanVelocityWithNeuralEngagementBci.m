function kalmanVelocityWithNeuralEngagementBci(controlCompSocket, expParams, okelecs)

global params codes

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

trialStarted = false;
timePtStarted = [];
timePtEnded = [];
bciStart = false;
timePtBciStarted = [];
timePtBciEnd = [];
samplesPerSecond = params.neuralRecordingSamplingFrequencyHz;%30000;
binSizeMs = expParams.binSizeMs;%50;
nasNetwork = expParams.nasNetwork;
velocity = [0; 0];
[nasNet.w1, nasNet.b1, nasNet.w2, nasNet.b2] = loadNasNet(nasNetwork);
gamma = expParams.gamma;
numBinsPrevNeZscAvg = expParams.numBinsPrevNeZscAvg;
neuralEngagementValueZscoredPrev = zeros(1, numBinsPrevNeZscAvg);

decoderParameterLocation = params.bciDecoderBasePathBciComputer;

msPerS = 1000;
samplesPerBin = binSizeMs/msPerS*samplesPerSecond;

binSpikeCountOverall = zeros(length(okelecs),1);
binSpikeCountNextOverall = zeros(length(okelecs), 1);

delValues = '';
tmstpInit = [];
waveforms = [];

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
    a = tic;
    [count,tmstpPrlEvt,events]=xippmex('digin');
    prlEvents = [events.parallel];
    tstpTrlStart = find(prlEvents==digitalCodeTrialStart);
    tstpTrlEnd = find(prlEvents==digitalCodeTrialEnd);
    tstpBciStart = find(prlEvents==digitalCodeBciStart);
    if length(tstpBciStart) > 1
        keyboard
    end
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
        velocity = [0; 0];
        neuralEngagementValueZscoredPrev = zeros(1, numBinsPrevNeZscAvg);
    end
  
    if trialStarted
        if ~isempty(ctrlMsg)
            % check if there's a new decoder or if there's an adjusted
            % velocity! (Important for calibration steps!)
            if strcmp(ctrlMsg, 'decoderParameterFile')
                decoderParameterFileRelativePath = receiveMessageSendAck(controlCompSocket);
                decoderParameterFileRelativePath(decoderParameterFileRelativePath=='\') = '/';
%                 decoderParameterFileRelativePath = 'satchel/21-Mar-2022/KalmanBci_14-02-09.mat';
                decoderParameterFileFullPath = fullfile(decoderParameterLocation, decoderParameterFileRelativePath);
                modelParams = load(decoderParameterFileFullPath);
                M0 = modelParams.M0;
                M1 = modelParams.M1;
                M2 = modelParams.M2;
                interpolatedConditionMeans = modelParams.interpolatedConditionMeans;
                interpolatedNeuralEngagementAxes = modelParams.interpolatedNeuralEngagementAxes;
                interpolatedNeuralEngagementValueMeans = modelParams.interpolatedNeuralEngagementValueMeans;
                interpolatedNeuralEngagementValueStds = modelParams.interpolatedNeuralEngagementValueStds;
                d = modelParams.d;
                betaOrthNorm = modelParams.betaOrthNorm;
                zScoreSpikesMat = modelParams.zScoreSpikesMat;
                    
                if isfield(modelParams, 'rotMat')
                    rotMat = modelParams.rotMat;
                else
                    rotMat = eye(2);
                end
                goodChannelNums = modelParams.channelsKeep;
                goodChannelInds = ismember(okelecs, goodChannelNums);
                fprintf('loaded new parameters from %s\n', decoderParameterFileRelativePath)
            elseif strcmp(ctrlMsg, 'requestParameters')
                parameterName = receiveMessageSendAck(controlCompSocket);
                parameterValue = modelParams.(parameterName);
                sendMessageWaitAck(controlCompSocket, typecast(parameterValue,'uint8'));
            else
                velocity = typecast(uint8(ctrlMsg), 'double')';
%                 fprintf('received constrained velocity\n');
            end
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
                timePtBciStarted = timePtBciStarted(timePtBciStarted>timePtStarted);
                timePtBinStart = timePtBciStarted;
                binSpikeCountOverall = zeros(length(goodChannelNums), 1);
                binSpikeCountNextOverall = zeros(length(goodChannelNums), 1);
                bciStart = true;
                bciJustStarted = true; % important for grabbing any early spikes
            end
            if ~isempty(tstpBciEnd)
                timePtBciEnd = tmstpPrlEvt(tstpBciEnd);
                timePtBciEnd = timePtBciEnd(timePtBciEnd>timePtStarted);
            end
            
            if timePtBciEnd > timePtBciStarted
                if bciStart
                    disp('bci end in trial')
                end
                bciStart = false;
                velocity = [0; 0];
                neuralEngagementValueZscoredPrev = zeros(1, numBinsPrevNeZscAvg);
            end
            
            if bciStart
                % Note to self: would it work to set okelecs to
                % goodChannels in some way?
%                 [~,tmstpInit, waveforms, ~]=xippmex('spike',okelecs,zeros(1,length(okelecs)));
                
                
                % only use data from channels that were good for
                % calibration
%                 tmstpGoodChSpk = tmstpInit(goodChannelInds);
%                 waveforms = waveforms(goodChannelInds);
%                 
%                 % run waveforms through NAS net
%                 tmstpNasSpk = cellfun(@(wvForms, tms) tms(runNASNetContinuous(w1, b1, w2, b2, wvForms, gamma)), waveforms, tmstpGoodChSpk, 'uni', 0);
% %                 tmstp = cellfun(@(wvForms, tms) tms, waveforms, tmstp, 'uni', 0);
%                 spikesThisBinByChannel = cellfun(@(x) x>=timePtBinStart & x<timePtBinStart+samplesPerBin, tmstpNasSpk, 'uni', 0);
%                 %                 waveformsThisBinByChannel = cellfun(@(wvFrm, spksInBin) wvFrm(spksInBin, :), waveforms, spikesThisBinByChannel, 'uni', 0);
%                 spikesNextBinByChannel = cellfun(@(x) x>=timePtBinStart+samplesPerBin, tmstpNasSpk, 'uni', 0);
%                 %                 waveformsNextBinByChannel = cellfun(@(wvFrm, spksInBin) wvFrm(spksInBin, :), waveforms, spikesNextBinByChannel, 'uni', 0);
%                 countsPerChannelCell = cellfun(@(x) sum(x, 2), spikesThisBinByChannel, 'uni', 0);
%                 countsExistChannel = ~cellfun('isempty', countsPerChannelCell);
%                 countsPerChannel = zeros(length(countsPerChannelCell), 1);
%                 countsPerChannel(countsExistChannel) = [countsPerChannelCell{countsExistChannel}];
%                 binSpikeCountOverall = binSpikeCountOverall + countsPerChannel;
                if bciJustStarted
                    [countsPerChannel, countsPerChannelNextBin] = countBinnedSpikesPerChannel(prevTmstpInit, goodChannelInds, timePtBinStart, samplesPerBin, nasNet, prevWaveforms, gamma);
                    binSpikeCountOverall = binSpikeCountOverall + countsPerChannel;
                    % this'll likely just be zero most of the time, but on bin
                    % edges it'll be needed
                    binSpikeCountNextOverall = binSpikeCountNextOverall + countsPerChannelNextBin;
                    bciJustStarted = false;
                end                
                [countsPerChannel, countsPerChannelNextBin] = countBinnedSpikesPerChannel(tmstpInit, goodChannelInds, timePtBinStart, samplesPerBin, nasNet, waveforms, gamma);
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
                
                % by checking whether there are timestamps from the *next*
                % bin, we confirm that we've completed the current bin and
                % can send off info
                if any(allTmstps>(timePtBinStart+samplesPerBin))
                    if ~any(allTmstps<=(timePtBinStart+samplesPerBin))
                        fprintf('(samples, time) past bin end: (%d, %d ms)\n', min(allTmstps-(timePtBinStart+samplesPerBin)), max(allTmstps-(timePtBinStart+samplesPerBin))/samplesPerSecond*msPerS);
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
                    [countsPerChannel, countsPerChannelNextBin] = countBinnedSpikesPerChannel(tmstpInit, goodChannelInds, timePtBinStart, samplesPerBin, nasNet, waveforms, gamma);
                    % scoop up last current bin spikes and grow the next
                    % bin spikes
                    binSpikeCountOverall = binSpikeCountOverall + countsPerChannel;
                    binSpikeCountNextOverall = binSpikeCountNextOverall + countsPerChannelNextBin;
                    
                    meanSpikeCount = mean(binSpikeCountOverall,2);
                    
%                     rrMat = [cosd(135) -sind(135); sind(135) cosd(135)];
%                     rotMat = rrMat * rotMat;
                    
                    
                    velocity = rotMat*(M0 + M1*rotMat'*velocity + M2 * meanSpikeCount);
                    
                    % Computing neural engagement value
                    
                    % Find angle of velocity (from 0 to 360)
                    intendedAngle = atan2d(velocity(2),velocity(1));
                    if((intendedAngle<0))
                        intendedAngle = intendedAngle + 360; % transforms domain from [-180,180] to [0,360]
                    end
                    intendedAngle = round(intendedAngle); % makes into integer
                    if(intendedAngle == 360)
                        intendedAngle = 0;
                    end
                    
                    % Choose the corresponding neural engagement axis
                    neuralEngagementAxisForCurrAng  = interpolatedNeuralEngagementAxes(:,intendedAngle+1);
                    conditionMeanForCurrAng         = interpolatedConditionMeans(:,intendedAngle+1);
                    
                    % Project spikes into orthonormalized FA space
                    factorActivity = betaOrthNorm * ((zScoreSpikesMat * meanSpikeCount)  - d);
                    
                    % Subtract condition-specific mean activity and project
                    % onto the corresponding neural engagement axis
                    neuralEngagementValue           = (factorActivity - conditionMeanForCurrAng)'*neuralEngagementAxisForCurrAng;
                    neuralEngagementValueZscored    = (neuralEngagementValue - interpolatedNeuralEngagementValueMeans(intendedAngle+1))...
                                                        /(interpolatedNeuralEngagementValueStds(intendedAngle+1));
                    
                    neuralEngagementValueZscoredTimeAvg = mean([neuralEngagementValueZscoredPrev(1:numBinsPrevNeZscAvg), neuralEngagementValueZscored]);
                    
%                     neuralEngagementValueZscoredTimeAvg
%                     neuralEngagementValueZscoredPrev
                    
                    neuralEngagementValueZscoredPrev = [neuralEngagementValueZscored neuralEngagementValueZscoredPrev(1:numBinsPrevNeZscAvg-1)];
                    

                    velocityAndNeuralEngagement   = velocity;
                    velocityAndNeuralEngagement(3)   = neuralEngagementValueZscoredTimeAvg;
                    
                    uint8Msg = typecast(velocityAndNeuralEngagement, 'uint8');
                    if any(uint8Msg==0)
                        %                     disp(uint8Msg)
                        %                     disp(typecast(meanSpikeCount, 'int8'))
                    end
                    msgToSend = uint8Msg';
%                     disp(velocity')
                    
                    matlabUDP2('send',controlCompSocket.sender, msgToSend);
%                     fprintf('sent unconstrained velocity\n');
                    
                     % the current bin is now what was the next bin before
                    binSpikeCountOverall = binSpikeCountNextOverall;
                    % zero out the counts for the next bin
                    binSpikeCountNextOverall(:) = 0;
                    timePtBinStart = timePtBinStart+samplesPerBin;
                end
            end
        end
    end
    loopTmAll = toc(loopTm);
end