%% Adaptive Sampling Closed Loop Code
% think about how to get trial start and end timing correct so we don't
% have empty spike counts
% zero spike counts in bcimat
% account for slow drift
% plot pcs
% plot trajectories
% pepe decode

close all; clear; clc;
addpath ../
addpath(genpath('../ex_control/'))
addpath(genpath('../xippmex-1.11'))
%addpath(genpath('../xippmex-1.4.0-test1/'));
addpath(genpath('../ex_bci/'))
%% Flags
% xippmex headstage setting. Either 'micro' or 'nano'
microornano= 'nano';

% sorting flag (only one of these can be equal to 1 or you will get an error
sortflag = 0;
crossnoiseflag = 0;
deepasortflag = 1;

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
filepath = '../../../../bciData/pepelepew/';


timebin = 0.0005;
numevents = 10;
loopsperbin = 8;

%% file parameter file 
%filepath = 'calib_online/hmmCalibSaves/';
%filename = 'Pe180815_s405a_dirmem_bci_discrete_0001_lowdHmm';
filename = [];
if isempty(filename)
    prethismonkey = strsplit(filepath,'/');
    thismonkey = prethismonkey{end-1}(1:2); %folder has format D:monkeyname
    thismonkey(1) = upper(thismonkey(1));
    thisyear = num2str(year(date));
    thisyear = thisyear(3:4);
    thismonth = num2str(month(date),'%02d');
    thisday = num2str(day(date),'%02d');
    thefile = dir([filepath,thismonkey,thisyear,thismonth,thisday,'*.mat']);
    
    if length(thefile)>1
        for n = 1:length(thefile)
            [~,thisfilename,~] = fileparts(thefile(n).name);
            userChoice = nan;
            while ~(strcmpi(userChoice,'y') || strcmpi(userChoice,'n'))
                userChoice = input([thisfilename,'.mat\nThis file? [y/n] ', ],'s');
            end
            
            switch lower(userChoice)
                case 'y'
                    filename = thisfilename;
                    break
                case 'n'
                otherwise
                    error('This shouldnt ever throw an error...something weird is going on');
            end
        end
        if isempty(filename)
            error('Please rerun code with desired filename as first input (do not include .nev extension)')
        end
    else
        [~,filename,~] = fileparts(thefile.name);
        fprintf(['Using the following file: \n',filename,'.mat\n'])  
        fprintf('If you want a different file, enter the filename\n as the first input to this function (do not include .nev extention)\n')
    end
end
thisyear = year(date);
thismonth = month(date);
thisday = day(date);
if str2double(['20' filename(3:4)])~=thisyear || str2double(filename(5:6))~=thismonth ||str2double(filename(7:8))~=thisday
    error('Wrong decoder parameter file!!!')
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Put Decoder/Plotting Variables Here: %%%%%%%%%%%%%%%%%%%%%
load([filepath filename]) % rename this file
keepneurons = modelparams.goodneuron;
if recalibflag == 1
    modelparams.recalibflag = 1;
else
    modelparams.recalibflag = 0;
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
postP = [];
trialstart = 0;
trialstartflag = 0;
num_bins = 1;
clockTime = clock;
meanvar1.mu = 0;
meanvar2.mu = 0;
meanvar1.sig = 1;
meanvar2.sig = 1;
anglein = 0;
distance = 0;
nCueTrials = 5;
nCorrectCueTrials = 0;
cueTrialCountBuffer = [];
correcttrials = [];
correctflag = 0;
numtrials = 1;
correcttrials = nan(10000,1);
percentcorrect = 0;
testcursor=0;
outputTimes = [];
outputTimesInd = 1;
aligntime = -1;
trialnum = 0;
firstspike = -1*ones(1,length(elecs));

%% parameters for deepasort
% constants
decisionboundary = 0.2;

% load parameter files
if deepasortflag==1
    sortdirectory = '../../../smithlab/matlab/neural_net_sort_stuff/';
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
% subplot(1,2,1)
% hold on
% mu = modelparams.meanproj;
% sig = modelparams.varproj;
% colors = distinguishable_colors(size(mu,2));
% for n = 1:size(mu,2)
%     func_plotEllipse(mu(:,n),sig{n},colors(n,:));
% end
% hproj=scatter(0,0,50,'r','filled')
% hold off
[~,tempname]=fileparts(filename);
%tempname = 'test80218';%
postind = 1;
testfile = ['../../../../bciData/pepelepew/',tempname,'_log',num2str(postind),'.txt'];
while exist(testfile,'file')
    postind = postind+1;
    testfile = ['../../../../bciData/pepelepew/',tempname,'_log',num2str(postind),'.txt'];
end
fileID = fopen(testfile,'a');
fprintf(fileID,'%s\t%s\t%s\t%s\t%s','trialnum','trialstarttime','binstarttimes', 'binendtimes', 'numspikesinbin');

if savespiketimesflag == 1
   fprintf(fileID,'\t%s','spiketimes');
end 

fprintf(fileID,'\t%s','sortedcount');

if recalibflag == 1
    fprintf(fileID,'\t%s','centeringmean');
end
fprintf(fileID,'\n');
fclose(fileID);

if videoflag == 1
    v = VideoWriter([filename,'bcimovie.avi'],'MPEG-4');
    v.FrameRate = 10;
    v.Quality = 20;
    open(v);
end
countstime = GetSecs;
which xippmex
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
            
            countbuffer = [countbuffer goodcounts(keepneurons==1)];
            
            % print the counts after filtering, sorting, etc.
            fprintf(fileID,'%i',sum(goodcounts));
            %% The following loop prints the recalibration means (We may not need this if the deepa waveform classification works)
            if recalibflag == 1
                for n = 1:length(modelparams.centeringMean)
                    fprintf(fileID,'%.3f,',modelparams.centeringMean(n));
                end
            end
            
            %                 for n = 1:length(matspiketimes)
            %                     fprintf(fileID,'%s,',num2str(matspiketimes(n)));
            %                 end
            
            
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
                countstime = double(GetSecs);
                countbuffer = [];
                postP = [];
                trialstart = 0;
                trialstartflag = 1;
                curlocx = [0];
                curlocy = [0];
                correctflag = 0;
                testcursor = 0;
                sendtime = countstime + 0.04;
                trialnum = num2str(receivecommand(11:end));
                
            elseif strcmp(receivecommand,'corrcue')
                nCorrectCueTrials = nCorrectCueTrials + 1;
                % add counts from correct cue trial to buffer
                cueTrialCountBuffer = [cueTrialCountBuffer countbuffer(:,1:20)];
                if nCorrectCueTrials == nCueTrials
                    % recenter mean
                    cueTrialMeans = mean(cueTrialCountBuffer,2);
                    cueTrialVars = var(cueTrialCountBuffer,[],2);
                    meanRates = cueTrialMeans .* (1/.05);
                    badChans = meanRates > 200 | meanRates < 1 | cueTrialVars==0;
                    numBadChans = sum(badChans);
                    if numBadChans>10
                        error('TOO MANY WONKY CHANNELS!!!!!');
                    else
                        fprintf('%d of %d channels are wonky.\n',numBadChans,length(badChans));
                    end
                    if recalibflag == 1
                        calibAngMean = modelparams.calib_angMeans(:,modelparams.targetAngles==anglein);
                        meanShift = cueTrialMeans - calibAngMean;
                        origCalibMean = mean(modelparams.calib_angMeans,2);
                        modelparams.centeringMean = origCalibMean + meanShift;
                    end
                    % clear buffer of previous correct cue trials
                    cueTrialCountBuffer = [];
                    % reset correct cue trial count to 0
                    nCorrectCueTrials = 0;
                end
            elseif strcmp('handshake',receivecommand)
                matlabUDP2('send',sockind,num2str(1));
            elseif strcmp('pretrialstart',receivecommand)
                if offlinemodeflag == 0; [~,~,~,units]= xippmex('spike',elecs,[zeros(1,length(elecs))]);end
            elseif strcmp('trialend',receivecommand)
                trialstartflag = 0;
                trialstart = 0;
                if correctflag == 1
                    correcttrials(numtrials) = 1;
                else
                    correcttrials(numtrials) = 0;
                end
                numtrials = numtrials+1;
                
            elseif strcmp(receivecommand(1:length('angle')),'angle')
                anglein = str2double(receivecommand((length('angle')+1):end));
                
            elseif strcmp(receivecommand(1:length('distance')),'distance')
                distance = str2double(receivecommand((length('distance')+1):end));
                
            elseif strcmp(receivecommand(1:length('ncuetrials')),'ncuetrials')
                nCueTrials = str2double(receivecommand((length('ncuetrials')+1):end));
                
            end
            
        end
        
        %% add code for sending information on each loop
        if trialstart == 1
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%% Add cursor decode Code: %%%%%%%%%%%%%%
            %[cursorpos1, cursorpos2] = decodeNBPCursor(countbuffer(:,end),modelparams); %define variables
            %             This code is for a discrete decoder followed by smoother
            %             [cursorpos1, cursorpos2,projnew] = decodeProjMVNCursor(countbuffer(:,end),modelparams);
            %             %[meanvarproj1,projlocx] = smoothCounts(projnew(1),modelparams,meanvarproj1,projlocx);
            %             %[meanvarproj2,projlocy] = smoothCounts(projnew(2),modelparams,meanvarproj2,projlocy);
            %             projlocx = projnew(1);
            %             projlocy = projnew(2);
            %             [meanvar1,curlocx] = smoothCounts(cursorpos1,modelparams,meanvar1,curlocx); % define variables
            %             [meanvar2,curlocy] = smoothCounts(cursorpos2,modelparams,meanvar2,curlocy); % define variables
            
            %             this code is for a continuous decoder (KF)
            %[meanvar1] = smoothCounts(countbuffer(:,end),modelparams.KFparamsdec,meanvar1);
            
            
            %%%% KF stuff %%%%
            %             if testcursor == 0
            %                 [meanvar] = ps9applyKF(countbuffer(:,end),modelparams.KFparamsdec);
            %             else
            %                 [meanvar] = ps9applyKF(countbuffer(:,end),modelparams.KFparamsdec,meanvar);
            %             end
            %             theta = deg2rad(anglein);
            %             targlocx = round(distance*cos(theta));
            %             targlocy = round(distance*sin(theta));
            %             curlocx = [curlocx meanvar.mu(1)];
            %             curlocy = [curlocy meanvar.mu(2)];
            %             sendcursorx = cumsum(curlocx);
            %             sendcursory = cumsum(curlocy);
            %             if offlinemodeflag==0
            %                 outpercent = (testcursor/10);
            %                 if double(GetSecs)-sendtime>bciperiod
            %                     fprintf('send loop slow %g\n',double(GetSecs)-sendtime)
            %                 end
            %                 while double(GetSecs)-sendtime<bciperiod
            %                 end
            %                 sendtime = countstime+0.04;
            %                 %tranmsit = 1
            %                 if demoflag == 1
            %                     matlabUDP2('send',sockind,num2str(round(outpercent*[targlocx targlocy])));
            %                 else
            %                     matlabUDP2('send',sockind,num2str(round([sendcursorx(end) sendcursory(end)])));
            %                 end
            %                 testcursor = testcursor + 1;
            %             end
            %%%%% end KF Stuff %%%%%
            
            %%%% start PNB stuff %%%%
            % [likelihoods] = likefuncPNB(neural) %actually grab posterior
            % sendposteriors
            if isempty(postP)
                %posterior = computePostNBP(sum(countbuffer(:,(end-loopsperbin+1):end),2),modelparams.mean);
                currPostP = hmm_onlineFilter(countbuffer(:,end),postP,modelparams);
                postP = [postP currPostP];
            else
                %posterior = computePostNBP(sum(countbuffer(:,:),2),modelparams.mean);
                currPostP = hmm_onlineFilter(countbuffer(:,end),postP(:,end),modelparams);
                postP = [postP currPostP];
            end
            if offlinemodeflag==0
                outpercent = (testcursor/10);
                if double(GetSecs)-sendtime>bciperiod
                    fprintf('send loop slow %g\n',double(GetSecs)-sendtime)
                end
                while double(GetSecs)-sendtime<bciperiod
                end
                sendtime = countstime+0.04;
                if demoflag == 1
                    targangleind = find(modelparams.allangles==anglein); %% need to add allangles to the model file
                    if targangleind <length( modelparams.allangles);sideind1 = targangleind+1;else sideind1 = 1;end
                    if targangleind >1;sideind2 = targangleind-1;else sideind2 = length(modelparams.allangles);end
                    demo = 0.1*ones(1,length(modelparams.allangles));
                    demo(targangleind) = outpercent;
                    demo(sideind1) = outpercent * 0.5 + 0.1*0.5;
                    demo(sideind2) = outpercent * 0.5 + 0.1*0.5;
                    if max(demo) > 1
                        demo = demo/max(demo);
                    end
                    matlabUDP2('send',sockind,num2str(demo,3));
                else
                    matlabUDP2('send',sockind,num2str(currPostP',3));
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
