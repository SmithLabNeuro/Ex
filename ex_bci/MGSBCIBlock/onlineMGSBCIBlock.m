%% Adaptive Sampling Closed Loop Code
% think about how to get trial start and end timing correct so we don't
% have empty spike counts
% zero spike counts in bcimat
% account for slow drift
% plot pcs
% plot trajectories
% pepe decode
%
close all; clear; clc;
addpath ../
addpath ../bciHMM
addpath bcihelpers
addpath(genpath('../../ex_control/'))
addpath(genpath('../../xippmex-1.11'))
%addpath(genpath('../xippmex-1.4.0-test1/'));
addpath('/home/smithlab/Dropbox/smithlab/matlab/fa')
%% make temp directory for online dat storage
bcitemploc = '/home/smithlab/Documents/bcitemp';%'../../../../../Documents/bcitemp';
mkdir(bcitemploc);

%% Flags


% xippmex headstage setting. Either 'micro' or 'nano'
microornano= 'nano';

% sorting flag (only one of these can be equal to 1 or you will get an error
sortflag = 0;
crossnoiseflag = 0;
deepasortflag = 0;

% recalibration
recalibflag = 0;

% data saving
savespiketimesflag = 0;

% multimedia flags
usesoundcardflag = 0;
videoflag = 0;

%simulation flags
offlinemodeflag = 0;
plotflag = 0;
demoflag = 0;

%% Initialize Loop Variables
bciperiod = 0.04996; %time in seconds
offbciperiod = 0.001;
looptimemaxerror = 0.01;
filepath = '/home/smithlab/bciData/satchel/';


timebin = 0.0005;
numevents = 10;
loopsperbin = 8;

%% file parameter file
%filepath = 'calib_online/hmmCalibSaves/';
%filename = 'Pe180815_s405a_dirmem_bci_discrete_0001_lowdHmm';
filename = [];%'Pe180920_s441a_distanceStabilityBCIcalib_0001_lowdHmm0'%
if isempty(filename)
    filename = findBCIFileName(filepath,'Sa');
end
thisyear = year(date);
thismonth = month(date);
thisday = day(date);
if str2double(['20' filename(3:4)])~=thisyear || str2double(filename(5:6))~=thismonth ||str2double(filename(7:8))~=thisday
    error('Wrong decoder parameter file!!!')
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Put Decoder/Plotting Variables Here: %%%%%%%%%%%%%%%%%%%%%
load([filepath filename]) % rename this file
%filename = 'testStability'

keepneurons = [];
modelparamsind = 1;
while isempty(keepneurons)
    keepneurons = modelparams(modelparamsind).goodneuron;% if two bcis, assume use same neurons. however it is possible that there may be empty modelparams, so need to find non-empty
    modelparamsind = modelparamsind + 1;
end
for n = 1:length(modelparams)
    if recalibflag == 1
        modelparams(n).recalibflag = 1;
    else
        modelparams(n).recalibflag = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize Xippmex and matlabUDP
if offlinemodeflag == 0
    matlabUDP2('close',0);
    sockind = matlabUDP2('open','192.168.2.10','192.168.2.11',4244);
    %addpath ../../smithlabrig/Ex/xippmex/
    status = xippmex;
    okelecs = xippmex('elec',microornano);
    elecs = okelecs;
    %matlabUDP2('send',sockind,num2str(1)); % send a wake up reminder if esc'd out
end

%% Initialize Psych Sounds
if usesoundcardflag == 1
    freq = 48000;
    reqlatency = 0;
    minLatency = 10 / 1000;
    AssertOpenGL;
    InitializePsychSound(1);
    PsychPortAudio('Verbosity', 0);
    paoutput = PsychPortAudio('Open', [], 1, 4, freq, 2, [], minLatency);
    outbuffersize = floor(freq*bciperiod);
    PsychPortAudio('FillBuffer', paoutput, zeros(2, outbuffersize));
    playbackstart = PsychPortAudio('Start', paoutput, 0, 0, 1);
    
    % initialize dummy sound variable
    soundcarddummy = zeros(1,length(1:bciperiod*freq));
    soundcarddummyoff = zeros(1,length(1:offbciperiod*freq));
    soundoncommand = 0;
end

%% Initialize keyboard
KbReleaseWait;
KbName('UnifyKeyNames');
RestrictKeysForKbCheck(KbName('ESCAPE'));

%% Initialize variables (need to clean this up)
if offlinemodeflag == 0;[count]=xippmex('spike',elecs,[zeros(1,length(elecs))]);end
countbuffer = [];
trialstart = 0;
trialstartflag = 0;
num_bins = 1;
x_new =[];
correctflag = -1;
numtrials = 1;
correcttrials = [];
percentcorrect = 0;
testcursor=0;
aligntime = -1;
trialnum = 0;
bcinum = 1; % this and all bci num instances will be needed if we want to do more than one decoder in a single session
postP = []; % this is the posterior probability vector for hmm smoothing.
%firstspike = -1*ones(1,length(elecs));
% for n = 1:length(modelparams)
%     if deltabciflag == 1
%         modelparams(n).distancemetric = @(x1,x2,mu1,mu2)(deltaBCIDistMetric(x1,x2,mu1,mu2));
%     else
%         modelparams(n).distancemetric = @(x,mu)((sum((x-repmat(mu,1,size(x,2))).^2,1))); % For some reason matlab won't let me use the original function in modelparams
%     end
% end
%% parameters for deepasort
% constants


% load parameter files
if deepasortflag==1
    decisionboundary = [];
    modelparamsind = 1;
    while isempty(decisionboundary)
        decisionboundary = modelparams(modelparamsind).gamma; % loop through to ensure that gamma field is non-empty. Assume all non-empty gamma are the same
        modelparamsind = modelparamsind + 1;
    end
    sortdirectory = '../../../../smithlab/matlab/neural_net_sort_stuff/';
    W1=dlmread([sortdirectory,'OldMonkey_Late_UberNet_w1'],' ');
    B1=dlmread([sortdirectory,'OldMonkey_Late_UberNet_b1'],' ');
    W2=dlmread([sortdirectory,'OldMonkey_Late_UberNet_w2'],' ');
    B2=dlmread([sortdirectory,'OldMonkey_Late_UberNet_b2'],' ');
    B1 = repmat(B1,1,1024);% the number of columns is set the buffersize of xippmex
    % find decision boundary
    sortboundary = log(decisionboundary/(1-decisionboundary))+B2(2)-B2(1);
end


%% init plot
f=figure;

%% init log file
[~,tempname]=fileparts(filename);
postind = 1;

testfile = [filepath,tempname,'_log',num2str(postind),'.txt'];
while exist(testfile,'file')
    postind = postind+1;
    testfile = [filepath,tempname,'_log',num2str(postind),'.txt'];
end
fileID = fopen(testfile,'a');
fprintf(fileID,'%s\t%s\t%s\t%s\t%s','trialnum','trialstarttime','binstarttimes', 'binendtimes', 'numspikesinbin');

if savespiketimesflag == 1
    fprintf(fileID,'\t%s','spiketimes');
end

fprintf(fileID,'\t%s','sortedcount');

if recalibflag == 1
    fprintf(fileID,'\t%s','meanparam');
end

fprintf(fileID,'\t%s','shamtrialind');

fprintf(fileID,'\n');
fclose(fileID);

%% init video if needed
if videoflag == 1
    v = VideoWriter([filename,'bcimovie.avi'],'MPEG-4');
    v.FrameRate = 10;
    v.Quality = 20;
    open(v);
end
countstime = GetSecs;
which xippmex % this is here so that we can confirm that we are using the correct xippmex version
%% BCI Loop
try
    while ~KbCheck
        % start loop time check
        if trialstart == 0
            countstime= GetSecs;
        end
        %% grab spike counts
        if trialstart == 1
            if offlinemodeflag == 0
                thistic = tic;
                [countbuffertemp,spiketimes,waveforms,units]=xippmex('spike',elecs,[zeros(1,length(elecs))]);
                %display(toc(thistic));
                %display(countbuffertemp(1))
                %                    minvals = cellfun(@(x)(min(x)),spiketimes,'UniformOutput',false);
                %                 mintime = min(cell2mat(minvals(cell2mat(cellfun(@(x)(~isempty(x)),minvals,'UniformOutput',false)))));
                %                 maxvals = cellfun(@(x)(max(x)),spiketimes,'UniformOutput',false);
                %                 maxtime = max(cell2mat(maxvals(cell2mat(cellfun(@(x)(~isempty(x)),maxvals,'UniformOutput',false)))));
                %
                % trialnum aligntime starttime endtime totalspikes
                %                 firstspikelist = nan(length(spiketimes),1);
                %                 for channelinds = 1:length(spiketimes)
                %                     if ~isempty(spiketimes{channelinds})
                %                         firstspikelist(channelinds) = spiketimes{channelinds}(1);
                %                     end
                %                 end
                % mintime = nanmin(firstspikelist(firstspikelist~=0))/30000;
                matspiketimes = [spiketimes{:}];
                mintime = min(matspiketimes(matspiketimes>0));
                if ~isempty(find(matspiketimes == 0,1))
                    warning('invalid spike times observed')
                end
                maxtime = max(matspiketimes);
                if mintime == 0
                    break
                end
                %maxtime = nanmax(firstspikelist(firstspikelist~=0))/30000;
                totalspikes = sum(countbuffertemp(:,end));
                fileID = fopen(testfile,'a');
                %fprintf(fileID,'%s\t%s\t%s\t%s\t%s\n',num2str(trialnum),num2str(aligntime),num2str(mintime), num2str(maxtime), num2str(totalspikes));
                fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t',num2str(trialnum),num2str(aligntime),num2str(mintime), num2str(maxtime), num2str(totalspikes));
                
                %% The following for loop prints min and max times for each channel (Do we need this? Need to check xippmex again)
                
                if savespiketimesflag == 1
                    for n = 1:length(elecs)
                        thisspikes = spiketimes{n};
                        if ~isempty(thisspikes)
                            thisspikegzero = thisspikes(thisspikes>0);
                            if length(thisspikegzero)<length(thisspikes)
                                warning('spikes with time equal 0 detected')
                            end
                            fprintf(fileID,'%i:%i,',(min(thisspikegzero)),(max(thisspikegzero)));
                        else
                            fprintf(fileID,'-1:-1,');
                        end
                    end
                    fprintf(fileID,'\t');
                end
                
            else
                countbuffertemp = spikingdata(num_bins).counts;
                spiketimes = spikingdata(num_bins).spiketimes;
                units = spikingdata(num_bins).units;
            end
            if sortflag ==1 && crossnoiseflag==0 && deepasortflag==0
                goodcounts = cellfun(@(x)(sum(x~=0)),units);
            elseif sortflag ==0 && crossnoiseflag==1 && deepasortflag==0
                [goodcounts, percentexclude] = excludeCrossChannelNoise(spiketimes,timebin,numevents);
            elseif sortflag ==0 && crossnoiseflag==0 && deepasortflag==1
                goodcounts = zeros(size(countbuffertemp));
                for n = 1:length(goodcounts)
                    thiswaves = round(waveforms{n}*10)/10; % this ensures that we have same format off and online
                    if ~isempty(thiswaves)
                        numwaves = size(thiswaves,1);
                        temp = max(0,thiswaves*W1+B1(:,1:numwaves)')*W2;
                        if numwaves == 1
                            goodcounts(n)=sum((temp(1)-temp(2))>sortboundary);%this approximates sum((W2'*max(0,W1*waveforms{n}'+b1))>4);%=
                        else
                            goodcounts(n)=sum((temp(:,1)-temp(:,2))>sortboundary);%this approximates sum((W2'*max(0,W1*waveforms{n}'+b1))>4);%
                        end
                    else
                        goodcounts(n)=0;
                    end
                end
            elseif sortflag ==0 && crossnoiseflag==0 && deepasortflag==0
                goodcounts = countbuffertemp;
            else
                error('Too many sort options selected');
            end
            
            countbuffer = [countbuffer goodcounts(modelparams(bcinum).goodneuron==1)];
            
            % print the counts after filtering, sorting, etc.
            fprintf(fileID,'%i',sum(goodcounts));
            %% The following loop prints the recalibration means (We may not need this if the deepa waveform classification works)
            %             if recalibflag == 1
            %                 fprintf(fileID,'\t');
            %                 for n = 1:length(modelparams(bcinum).mu)
            %                         fprintf(fileID,'%.3f,',modelparams(bcinum).mu(n));
            %                 end
            %             end
            
            %                 for n = 1:length(matspiketimes)
            %                     fprintf(fileID,'%s,',num2str(matspiketimes(n)));
            %                 end
            %   fileID = fopen(testfile,'a');
            fprintf(fileID,'\t');
            fprintf(fileID,'\n');
            fclose(fileID);
        end
        
        %% check for messages
        if offlinemodeflag == 1
            messages = spikingdata(num_bins).messages;
            messagetotal = length(messages);
            currmessagecount = 1;
            num_bins = num_bins + 1;
            if num_bins >length(spikingdata)
                break
            end
        end
        
        while ((offlinemodeflag==1)&&(currmessagecount<=messagetotal)) || ((offlinemodeflag==0)&&matlabUDP2('check',sockind))
            if offlinemodeflag == 0
                receivecommand = matlabUDP2('receive',sockind);
            else
                receivecommand = messages{currmessagecount};
                currmessagecount = currmessagecount + 1;
            end
            if length(receivecommand)>10&&strcmp('trialstart',receivecommand(1:10))
                if offlinemodeflag == 0
                    aligntime = xippmex('time');
                    [~,~,~,units]= xippmex('spike',elecs,[zeros(1,length(elecs))]);
                    
                end
                %fclose(fileID);
                countstime = double(GetSecs);
                countbuffer = [];
                trackvals = [];
                trialstart = 0;
                trialstartflag = 1;
                correctflag = -1;
                testcursor = 0;
                postP = [];
                sendtime = countstime + 0.04;
                trialnum = num2str(receivecommand(11:end));
            elseif strcmp('handshake',receivecommand)
                matlabUDP2('send',sockind,num2str(1));
                
            elseif length(receivecommand)>=11&&strcmp('recalibrate',receivecommand(1:11))
                %startRecalibTrial = str2double(receivecommand(12:end));
                %modelparamsold = modelparams(1); % only need to put in one of the bci params. Note: this line is a bad assumption: should not assume modelparams(1) is non-empty
                %[modelparams] = calibrateDistanceBCIFA(filepath,filepath,'filename',[filename,'_'],'dat',onlinedat,'modelparamsold',modelparamsold,'trialnum',numtrials,'dimsbci1',modelparams(1).dimsbci,'dimsbci2',modelparams(2).dimsbci,'deltabciflag',1,'targetCorrect',0.3); % need to fix the last two entries, these are bad assumptions
                %fprintf('Calibration finished!\n')
                fprintf('Could not recalibrate...\n')
            elseif strcmp('pretrialstart',receivecommand)
                if offlinemodeflag == 0; [~,~,~,units]= xippmex('spike',elecs,[zeros(1,length(elecs))]);end
            elseif length(receivecommand)>=8&&strcmp('trialend',receivecommand(1:8))
                if length(receivecommand)>8
                    correctflag = str2double(receivecommand(9))-1;
                end
                trialstartflag = 0;
                trialstart = 0;
                numcalibtrials = 20;
                if correctflag == 1
                    correcttrials(numtrials) = 1;
                    onlinedat(numtrials).counts = countbuffer;
                    onlinedat(numtrials).cursorloc = trackvals;
                    onlinedat(numtrials).correctflag = 1;
                    onlinedat(numtrials).bcinum = bcinum;
                    onlinedat(numtrials).trialnum = numtrials;
                    onlinedat(numtrials).postP = postP
                    temponlinedat = onlinedat(numtrials);
                    save([bcitemploc,'/onlinedat',num2str(numtrials)],'temponlinedat')
                    numtrials = numtrials+1;
                elseif correctflag == 0
                    correcttrials(numtrials) = 0;
                    onlinedat(numtrials).counts = countbuffer;
                    onlinedat(numtrials).cursorloc = trackvals;
                    onlinedat(numtrials).correctflag = 0;
                    onlinedat(numtrials).bcinum = bcinum;
                    onlinedat(numtrials).trialnum = numtrials;
                    onlinedat(numtrials).postP = postP;
                    temponlinedat = onlinedat(numtrials);
                    save([bcitemploc,'/onlinedat',num2str(numtrials)],'temponlinedat')
                    numtrials = numtrials+1;
                    modelparamsold = modelparams;
                    %                     if recalibflag==1&&length(correcttrials)>numcalibtrials && sum(correcttrials((end-numcalibtrials+1):end))==0
                    %                         modelparams = calibrateDistanceBCIFA([filename,'_'],filepath,filepath,modelparamsold.gamma,modelparamsold.targetCorrect,1,modelparams.multgain,numtrials,onlinedat((end-numcalibtrials+1):end),modelparamsold);
                    %                         fprintf('Too Many Misses. Updating Model Params.\n')
                    %                     end
                end
                if exist('onlinedat','var')&&length(onlinedat)>60&&mod(length(onlinedat),10)==0&&(correctflag>=0)
                    thisonlinedat = onlinedat;
                    percorrectvec = mean([thisonlinedat.correctflag]);
                    fprintf('Overall Percent Correct: %g\n',percorrectvec)
                    if length(modelparams)>1
                        onlinedat1 = thisonlinedat([thisonlinedat.bcinum]==1);
                        if ~isempty(onlinedat1)
                            fprintf('BCI 1 Percent Correct: %g\n',mean([onlinedat1.correctflag]))
                        end
                        onlinedat2 = thisonlinedat([thisonlinedat.bcinum]==2);
                        if ~isempty(onlinedat2)
                            fprintf('BCI 2 Percent Correct: %g\n',mean([onlinedat2.correctflag]))
                        end
                    end
                end
            end
        end
        
        %% add code for sending information on each loop
        if trialstart == 1
            
            %             this code is for a continuous decoder
            
            
            % new stuff
            if isempty(postP)
                %posterior = computePostNBP(sum(countbuffer(:,(end-loopsperbin+1):end),2),modelparams.mean);
                currPostP = hmm_onlineFilter(countbuffer(:,end),postP,modelparams);
                postP = [postP currPostP];
            else
                %posterior = computePostNBP(sum(countbuffer(:,:),2),modelparams.mean);
                currPostP = hmm_onlineFilter(countbuffer(:,end),postP(:,end),modelparams);
                postP = [postP currPostP];
            end
            
            
            [cursorX,cursorY] = posterior2pos(modelparams.allangles,postP(:,end));
            valtosend = [cursorX;cursorY]';
            
            trackvals = [trackvals valtosend];
            if offlinemodeflag==0
                outpercent = (testcursor/10);
                if double(GetSecs)-sendtime>bciperiod
                    fprintf('send loop slow %g\n',double(GetSecs)-sendtime)
                end
                while double(GetSecs)-sendtime<bciperiod
                end
                sendtime = countstime+0.04;
                if demoflag == 1
                    demo = outpercent;
                    matlabUDP2('send',sockind,num2str(demo,3));
                else
                    matlabUDP2('send',sockind,num2str(valtosend,20));
                end
                testcursor = testcursor + 1;
            end
            
            
            
            %%%% end PNB stuff %%%%
            
            
            %%%%%% Add Plot Code %%%%%%
            %             subplot(1,2,1)
            %             set(hproj,'XData',projlocx,'YData',projlocy)
            %             drawnow
            
            %             subplot(1,2,2)
            if plotflag == 1
                plot(sendcursorx,sendcursory,'Color','b')
                
                hold on
                scatter(targlocx,targlocy,10,'k','fill')
                scatter(0,0,1,'k')
                t = linspace(0,2*pi);plot(40*cos(t)+targlocx,40*sin(t)+targlocy)
                xlim([-140 140])
                ylim([-140 140])
                text(50,-100,['% Correct: ',num2str(100*percentcorrect,3)])
                hold off
                drawnow
            end
            if videoflag == 1
                frame = getframe(f);
                writeVideo(v,frame.cdata)
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        %% load sound to soundcard buffer
        if trialstartflag == 1
            trialstart = 1;
        end
        if offlinemodeflag == 0
            if trialstart == 1
                if usesoundcardflag == 1;audiodata = [soundcarddummy;soundcarddummy];end
                currlooptime = bciperiod;
            else
                if usesoundcardflag == 1;audiodata = [soundcarddummyoff;soundcarddummyoff];end
                currlooptime = offbciperiod;
            end
            if usesoundcardflag == 1;PsychPortAudio('FillBuffer', paoutput, audiodata, 1);end
        else
            currlooptime = offbciperiod;
        end
        %% check that timing is correct
        
        time = double(GetSecs) - countstime;
        % check if loop time is too long
        if trialstart==1 && (time-currlooptime)>looptimemaxerror*bciperiod
            display(testcursor)
            fprintf('Loop time is too long: %d\n',time)
            countstime = double(GetSecs);
        else
            % otherwise wait until loop time is finished
            
            while double(GetSecs)-countstime < currlooptime
            end
            %            if trialstart == 1
            %                tcounts = toc(countstime)
            %                tsend = toc(sendtime)
            %                difftime = toc(countstime) - toc(sendtime)
            %            end
            countstime = double(GetSecs);
        end
        %time = toc;
        %display(time)
        
        %% initialize next spikecount vec
        
        
    end
    
    %% clean up
    testfile = [filepath,filename,'onlinedat1.mat'];
    postind = 1;
    while exist(testfile,'file')
        postind = postind+1;
        testfile = [filepath,filename,'onlinedat',num2str(postind),'.mat'];
    end
    save(testfile,'onlinedat')
    
    if ~isempty(strfind(bcitemploc,'bcitemp')) % check that you don't accidentally delete files in the wrong folder
        delete([bcitemploc,'/onlinedat*.mat']);
    end
    
    
    
    if usesoundcardflag == 1
        PsychPortAudio('Stop', paoutput, 1);
        PsychPortAudio('Close');
        
    end
    if videoflag == 1
        close(v)
    end
    if offlinemodeflag == 0
        xippmex('close');
        matlabUDP2('all_close');
    end
    %% (Optional) save output
catch err
    %% clean up
    display('ERROR!!!')
    if usesoundcardflag == 1
        PsychPortAudio('Stop', paoutput, 1);
        PsychPortAudio('Close');
    end
    if offlinemodeflag == 0
        xippmex('close');
        matlabUDP2('all_close');
    end
    %% (Optional) save output
end
