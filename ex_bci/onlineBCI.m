%% Adaptive Sampling Closed Loop Code
% think about how to get trial start and end timing correct so we don't
% have empty spike counts
% zero spike counts in bcimat
% account for slow drift
% plot pcs
% plot trajectories
% pepe decode

clear all;
close all;
addpath(genpath('../'))
%% Flags
microornano= 'nano';
sortflag = 0;
crossnoiseflag = 0;
usesoundcardflag = 0;
offlinemodeflag = 0;
videoflag = 0;
plotflag = 0;
demoflag = 0;
%% Initialize Loop Variables
bciperiod = 0.04996; %time in seconds
offbciperiod = 0.001;
looptimemaxerror = 0.01;
keepneurons = ones(96,1);
filename = 'calib_online/Wa171201_s187a_dirmem_bci_0002modelparamstestgauss50msalltrialonlinetest.mat';
timebin = 0.0005;
numevents = 20;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Put Decoder/Plotting Variables Here: %%%%%%%%%%%%%%%%%%%%%
load([filename]) % rename this file
%modelparams = pmdmodelparams;
if offlinemodeflag == 1; load([filename,'spikingdata']); end %spikingdata = trainingdata; 
keepneurons = modelparams.goodneuron;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize Xippmex and matlabUDP
if offlinemodeflag == 0
    matlabUDP2('close',0);
    sockind = matlabUDP2('open','192.168.2.10','192.168.2.11',4244);
    addpath ../../smithlabrig/Ex/xippmex/
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
clockTime = clock;
meanvar1.mu = 0;
meanvar2.mu = 0;
meanvar1.sig = 1;
meanvar2.sig = 1;
anglein = 0;
distance = 0;
correcttrials = [];
correctflag = 0;
numtrials = 1;
correcttrials = nan(10000,1);
percentcorrect = 0;
testcursor=0;
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

if videoflag == 1 
v = VideoWriter([filename,'bcimovie.avi'],'MPEG-4');
v.FrameRate = 10;
v.Quality = 20;
open(v);
end
countstime = GetSecs;
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
                [countbuffertemp,spiketimes,~,units]=xippmex('spike',elecs,[zeros(1,length(elecs))]);
            else
                countbuffertemp = spikingdata(num_bins).counts;
                spiketimes = spikingdata(num_bins).spiketimes;
                units = spikingdata(num_bins).units;
            end
            if sortflag ==1
                goodcounts = cellfun(@(x)(sum(x~=0)),units);
            elseif crossnoiseflag == 1
                [goodcounts, percentexclude] = excludeCrossChannelNoise(spiketimes,timebin,numevents);
            else
                goodcounts = countbuffertemp;
            end     
            countbuffer = [countbuffer goodcounts(keepneurons==1)];
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
            if strcmp('trialstart',receivecommand)
                if offlinemodeflag == 0; [~,~,~,units]= xippmex('spike',elecs,[zeros(1,length(elecs))]);end
                countstime = double(GetSecs);
                countbuffer = [];
                trialstart = 0;
                trialstartflag = 1;
                curlocx = [0];
                curlocy = [0];
                correctflag = 0;
                testcursor = 0;
                sendtime = countstime + 0.04;
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
                if numtrials>25
                    percentcorrect = nanmean(correcttrials(numtrials-25:numtrials));
                end
                fprintf('Percent Correct = %g\n',percentcorrect);
                numtrials = numtrials+1;
            elseif strcmp(receivecommand(1:length('angle')),'angle')
                anglein = str2double(receivecommand((length('angle')+1):end));
                
            elseif strcmp(receivecommand(1:length('distance')),'distance')
                distance = str2double(receivecommand((length('distance')+1):end));       
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
            if testcursor == 0
                [meanvar] = ps9applyKF(countbuffer(:,end),modelparams.KFparamsdec);
            else
                [meanvar] = ps9applyKF(countbuffer(:,end),modelparams.KFparamsdec,meanvar);
            end
            theta = deg2rad(anglein);
            targlocx = round(distance*cos(theta));
            targlocy = round(distance*sin(theta));
            curlocx = [curlocx meanvar.mu(1)];
            curlocy = [curlocy meanvar.mu(2)];
            sendcursorx = cumsum(curlocx);
            sendcursory = cumsum(curlocy);
            if offlinemodeflag==0
                outpercent = (testcursor/10);
                if double(GetSecs)-sendtime>bciperiod
                    fprintf('send loop slow %g\n',double(GetSecs)-sendtime)
                end
                while double(GetSecs)-sendtime<bciperiod
                end
                sendtime = countstime+0.04;
                %tranmsit = 1
                if demoflag == 1
                    matlabUDP2('send',sockind,num2str(round(outpercent*[targlocx targlocy])));
                else
                    matlabUDP2('send',sockind,num2str(round([sendcursorx(end) sendcursory(end)])));
                end
                testcursor = testcursor + 1;
            end
            
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
            if norm([min(abs(targlocx-sendcursorx)) min(abs(targlocy-sendcursory))],2)<40
                correctflag = 1;
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
