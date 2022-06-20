function photoDiodeSetBci(controlCompSocket, expParams, okelecs)


digitalCodeTrialStart = 1;
digitalCodeTrialEnd = 255;
tstpTrlStart = [];
trialStarted = false;
samplesPerSecond = 30000;
binSizeMs = 50;
msPerS = 1000;
samplesPerBin = binSizeMs/msPerS*samplesPerSecond;

spikeCountOverall = zeros(length(okelecs),1);
delValues = '';

%     pause(0.1)
[count,tmstp,events]=xippmex('digin');
prlEvents = [events.parallel];
disp(prlEvents);
% check the start of the trial
if ~trialStarted
    tstpTrlStart = find(prlEvents==digitalCodeTrialStart);
    if length(tstpTrlStart)>1
        disp('missed a trial')
        tstpTrlStart = tstpTrlStart(end);
    end
    if ~isempty(tstpTrlStart)
        disp('trial start')
        timePtStarted = tmstp(tstpTrlStart);
        timePtLastBin = timePtStarted;
        trialStarted = true;
        delValues = '';
    else
        % check the connection between trials...
        bciEnd = checkIfBciEnd(controlCompSocket);
        if bciEnd
            return
        end
    end
end

% accumulate spike count if trial has started
if trialStarted
    tstpTrlEnd = find(prlEvents==digitalCodeTrialEnd);
    if length(tstpTrlEnd)>1
        disp('missed a trial')
        tstpTrlEnd = tstpTrlEnd(end);
    end
    if ~isempty(tstpTrlEnd)
        trialStarted = false;
        timePtEnded = tmstp(tstpTrlEnd);
        [count,tmstp, waveforms, units]=xippmex('spike',okelecs,[zeros(1,length(okelecs))]);
        waveformsUse = cellfun(@(x) x>timePtStarted & x<timePtEnded, tmstp, 'uni', 0);
        countsPerChannel = cellfun(@(x) sum(x, 2), waveformsUse);
        spikeCountOverall = spikeCountOverall + countsPerChannel;
        strHere = sprintf('%4d ', spikeCountOverall');
        %             fprintf('%s%s', delValues, strHere);
        fprintf('\n')
        fprintf('mean trial spike count : %d\n', mean(spikeCountOverall));
        trialTimeInSeconds = (timePtEnded - timePtStarted)/samplesPerSecond;
        fprintf('mean FR (count/s) : %d\n', mean(spikeCountOverall)/(binSizeMs/msPerS));
        disp('trial end')
        delValues = '';
        spikeCountOverall = zeros(length(okelecs),1);
        %             trialSpikeCountAverage = [];
    else
        [count,tmstp, waveforms, units]=xippmex('spike',okelecs,[zeros(1,length(okelecs))]);
        waveformsUse = cellfun(@(x) x>timePtStarted, tmstp, 'uni', 0);
        countsPerChannel = cellfun(@(x) sum(x, 2), waveformsUse);
        
        allTmstps = cat(2, tmstp{:});
        if any(allTmstps>(timePtLastBin+samplesPerBin))
            meanSpikeCount = mean(spikeCountOverall);
            if isnan(meanSpikeCount)
                keyboard
            else
                disp(meanSpikeCount);
            end
            uint8Msg = typecast(meanSpikeCount, 'uint8');
            if any(uint8Msg==0)
                %                     disp(uint8Msg)
                %                     disp(typecast(meanSpikeCount, 'int8'))
            end
            msgToSend = uint8Msg;
            
            matlabUDP2('send',controlCompSocket, msgToSend);
            spikeCountOverall = zeros(length(okelecs),1);
            timePtLastBin = timePtLastBin+samplesPerBin;
        else
            spikeCountOverall = spikeCountOverall + countsPerChannel;
            strHere = sprintf('%4d ', spikeCountOverall');
            %                 fprintf('%s%s', delValues, strHere);
            
            delValues = repmat(sprintf('\b'), 1, length(strHere));
        end
        
        
    end
end