function bcistruct = createBCIStruct(filepath,hivepath,trainingfile,bcifile,bcilog,paramsfile)
% example:
% filestem = 'Wa180222_s234a_dirmem_bci_discrete_';
% bcistruct = createBCIStruct([],[filestem,'0001'],[filestem,'0002'],[filestem,'0001_lowdHmm_log1.txt'],[filestem,'0001_lowdHmm'])

% set filepath = [] if the files are in the current directory, else specify
% the directory.
%% set path
addpath(genpath('nev2bcistruct'))

%% constants
ALIGN = 253;
SOUND_CHANGE = 134;
%% load model params
disp('loading model parameters')
filecount = 1;
%if iscell(paramsfile)    
    for n = 1:length(paramsfile)
        %modelparamsstruct(n) = load([hivepath,paramsfile{n}]);
        thisfilestruct = dir([hivepath,paramsfile{n},'*.mat']);
        for m = 1:length(thisfilestruct)
            if ~contains(thisfilestruct(m).name,'onlinedat')
            if filecount == 1
                modelparamsstruct = load([hivepath,thisfilestruct(m).name]);
            else
                modelparamsstruct(filecount) = load([hivepath,thisfilestruct(m).name]);
            end
            filecount = filecount + 1;
            end
        end
    end
    % assume that all modelparams use the same gamma
    try
        gamma = modelparamsstruct(1).modelparams.gamma;
    catch
        error('model param file was empty: check that name of .nev file and input file name are correct')
    end
%else
%     modelparamsstruct = load([hivepath,paramsfile]);
%end

%% create training struct
disp('parsing training file')
if iscell(trainingfile)
    for n = 1:length(trainingfile)
        if n == 1
            %trainingdat = nev2bcistruct([filepath,trainingfile{n},'.nev']);
            [trainingdat,~] = nev2sortedStruct([filepath,trainingfile{n},'.nev'],gamma,false);
        else
            %trainingdat =[trainingdat nev2bcistruct([filepath,trainingfile{n},'.nev'])];
            [trainingdattemp, ~] = nev2sortedStruct([filepath,trainingfile{n},'.nev'],gamma,false);
            trainingdat = [trainingdat trainingdattemp];
        end
    end
else
    [trainingdat,nev] = nev2sortedStruct([filepath,trainingfile,'.nev'],gamma,false);
   %     nev = readNEV([filepath,trainingfile,'.nev']);
    calibrationnevstruct.nevchan = uint16(nev(nev(:,1)~=0,1));
    calibrationnevstruct.nevsort = uint16(nev(nev(:,1)~=0,2));
    calibrationnevstruct.firstspike = uint32(nev(1,3)*30000);
    calibrationnevstruct.timediff = uint16(diff(nev(nev(:,1)~=0,3)*30000));
    %calibrationnevstruct.time = uint32(nev(nev(:,1)~=0,3)*30000);
    calibrationnevstruct.nev_info = NEV_displayheader([filepath,trainingfile,'.nev']);
    for ii = 1:length(trainingdat)
        trainingdat(ii).nevinfo.nevclockstart = calibrationnevstruct.nev_info.nevclockstart;
    end

    
    %trainingdat = nev2bcistruct(nev,'nevreadflag',1,'nev_info',calibrationnevstruct.nev_info); 

    
end

%% create bci struct
disp('parsing bci file')
numbcitrials = [];
bcistructstarts = [];
bcistructends = [];
if iscell(bcifile)
    for n = 1:length(bcifile)
        if n == 1
            
            [bcidat,nev] = nev2sortedStruct([filepath,bcifile{n},'.nev'],gamma,false);
            %nev = readNEV([filepath,bcifile{n},'.nev']);
            bcinevstruct(n).nevchan = uint16(nev(nev(:,1)~=0,1));
            bcinevstruct(n).firstspike = uint32(nev(1,3)*30000);
            bcinevstruct.nevsort = uint16(nev(nev(:,1)~=0,2));
            bcinevstruct(n).timediff = uint16(diff(nev(nev(:,1)~=0,3)*30000));
            bcinevstruct(n).nev_info = NEV_displayheader([filepath,bcifile{n},'.nev']);
            %bcidat = nev2bcistruct(nev,'nevreadflag',1,'nev_info',bcinevstruct(n).nev_info); 
            %bcidat = nev2bcistruct([filepath,bcifile{n},'.nev']);
            for bciind = 1:length(bcidat)
                bcidat(bciind).bcitrial = bcidat(bciind).params.trial.bciTrial;
                bcidat(bciind).bciStarted = ~isempty(find(bcidat(bciind).trialcodes(:,2)==ALIGN,1));
                bcidat(bciind).nevinfo.nevclockstart = bcinevstruct(n).nev_info.nevclockstart;
            end
            bcitrials = [bcidat.bciStarted];
            numbcitrials = sum(bcitrials);
            bcistructstarts = 1;
            bcistructends = length(bcidat);
        else
            bcistructstarts = [bcistructstarts length(bcidat)+1];
            [tempbci,nev] = nev2sortedStruct([filepath,bcifile{n},'.nev'],gamma,false);
            %nev = readNEV([filepath,bcifile{n},'.nev']);
            bcinevstruct(n).nevchan = uint16(nev(nev(:,1)~=0,1));
            bcinevstruct.nevsort = uint16(nev(nev(:,1)~=0,2));
            bcinevstruct(n).firstspike = uint32(nev(1,3)*30000);
            bcinevstruct(n).timediff = uint16(diff(nev(nev(:,1)~=0,3)*30000));
            bcinevstruct(n).nev_info = NEV_displayheader([filepath,bcifile{n},'.nev']);
            %tempbci = nev2bcistruct(nev,'nevreadflag',1,'nev_info',bcinevstruct(n).nev_info); 
            %tempbci = nev2bcistruct([filepath,bcifile{n},'.nev']);
            for bciind = 1:length(tempbci)
                tempbci(bciind).bcitrial = tempbci(bciind).params.trial.bciTrial;
                tempbci(bciind).bciStarted = ~isempty(find(tempbci(bciind).trialcodes(:,2)==ALIGN,1));
                tempbci(bciind).nevinfo.nevclockstart = bcinevstruct(n).nev_info.nevclockstart;
            end
            bcitrials = [tempbci.bciStarted]; 
            numbcitrials = [numbcitrials sum(bcitrials)];
            bcidat =[bcidat tempbci];
            bcistructends = [bcistructends length(bcidat)];
        end
        
    end
else
    [bcidat,nev] = nev2sortedStruct([filepath,bcifile,'.nev'],gamma,false);
    %nev = readNEV([filepath,bcifile,'.nev']);
    bcinevstruct.nevchan = uint16(nev(nev(:,1)~=0,1));
    bcinevstruct.nevsort = uint16(nev(nev(:,1)~=0,2));
    bcinevstruct.firstspike = uint32(nev(1,3)*30000);
    bcinevstruct.timediff = uint16(diff(nev(nev(:,1)~=0,3)*30000));
    bcinevstruct.nev_info = NEV_displayheader([filepath,bcifile,'.nev']);
    %bcidat = nev2bcistruct(nev,'nevreadflag',1,'nev_info',bcinevstruct.nev_info); 
    for bciind = 1:length(bcidat)
        bcidat(bciind).bcitrial = bcidat(bciind).params.trial.bciTrial;
        bcidat(bciind).bciStarted = ~isempty(find(bcidat(bciind).trialcodes(:,2)==ALIGN,1));
        bcidat(bciind).nevinfo.nevclockstart = bcinevstruct.nev_info.nevclockstart;
    end
    bcitrials = [bcidat.bciStarted];
    numbcitrials = sum(bcitrials);
    bcistructstarts = 1;
    bcistructends = length(bcidat);
end

%% load bci log information
disp('loading log file')
if iscell(bcilog)    
    for n = 1:length(bcilog)
        if n == 1
            bcitimes = tdfread([hivepath,bcilog{n}]);
            logfields = fieldnames(bcitimes);
        else
            bcitimes2 =tdfread([hivepath,bcilog{n}]);
            for fieldind = 1:length(logfields)
                tempfield = logfields{fieldind};
                bcitimes.(tempfield) = [bcitimes.(tempfield);bcitimes2.(tempfield)];
            end
        end
    end
else
     bcitimes = tdfread([hivepath,bcilog]);
end


%% make 1 ms spike counts bins (entered in as field "counts_1ms")
disp('Making 1ms spike count bins')

%[~, ~, ~,dummy,~,~] = prepCalibCounts(trainingdat, 100, 3,0,0.001,10000);
%for n = 1:length(dummy); trainingdat(n).counts_1ms = sparse(logical(dummy(n).counts));end
%[~, ~, ~,dummy,~,~] = prepCalibCounts(bcidat, 100, 3,0,0.001,10000);
%for n = 1:length(dummy); bcidat(n).counts_1ms = sparse(logical(dummy(n).counts));end


%% find which series of log file lines match up with the bci structs
disp('aligning log file with trials')
startlines = [1; find(diff(bcitimes.trialnum)<0)+1];
if startlines(end) == length(bcitimes.trialnum)
    startlines(end) = [];
end
endlines = startlines(startlines>1)-1;
if bcitimes.trialnum(end) ~= 1
    endlines = [endlines; length(bcitimes.trialnum)];
end

goodstarts = zeros(length(startlines),1);
for n = 1:length(startlines)
    thislog = unique(bcitimes.trialnum(startlines(n):endlines(n)));
    for m = 1:length(numbcitrials)
        if  numbcitrials(m) == length(thislog)
            goodstarts(n) = m;
        end
    end
end

%if there is only one  log file and bci file
if length(startlines)==1&&length(bcistructstarts)==1 
    goodstarts = 1;
end

%% align log file and bci trial data
bcidat(1).trialstarttime = [];
bcidat(1).binstarttimes = [];
bcidat(1).binendtimes = [];
bcidat(1).numspikesinbin = [];
bcidat(1).centeringmean = [];
bcidat(1).sortedcount = [];
bcidat(1).shamtrialind = [];
for n = 1:length(bcistructstarts)
    n
    goodstartline = startlines(goodstarts==n);
    goodendline = endlines(goodstarts==n);
    fullthislog = bcitimes.trialnum(goodstartline:goodendline);
   % fullthislog(fullthislog >=1180) = fullthislog(fullthislog>=1180)-1;
    thistrialstarttime = bcitimes.trialstarttime(goodstartline:goodendline);
    thisbinstarttimes = bcitimes.binstarttimes(goodstartline:goodendline);
    thisbinsendtimes = bcitimes.binendtimes(goodstartline:goodendline);
    thisnumspikesinbin = bcitimes.numspikesinbin(goodstartline:goodendline);
    thissortedcount = bcitimes.sortedcount(goodstartline:goodendline);
    if isfield(bcitimes,'shamtrialind')
        thisshamtrialind = bcitimes.shamtrialind(goodstartline:goodendline);
    end
    if isfield(bcitimes,'centeringmean')
        this_centeringmean = bcitimes.centeringmean(goodstartline:goodendline,:);
    end
    thislog = unique(fullthislog);
    tempbci = bcidat(bcistructstarts(n):bcistructends(n));
    for m = 1:length(thislog)   
        indices = fullthislog==thislog(m);
        tempbci(thislog(m)).trialstarttime = thistrialstarttime(find(indices,1));
        tempbci(thislog(m)).binstarttimes = thisbinstarttimes(indices);
        tempbci(thislog(m)).binendtimes = thisbinsendtimes(indices);
        tempbci(thislog(m)).numspikesinbin = thisnumspikesinbin(indices);
        tempbci(thislog(m)).sortedcount = thissortedcount(indices);
        if isfield(bcitimes,'shamtrialind')
            tempbci(thislog(m)).shamtrialind = thisshamtrialind(indices);
        end
        if isfield(bcitimes,'centeringmean')
            tempbci(thislog(m)).centeringmean = this_centeringmean(find(indices,1),:);
        end
    end
    bcidat(bcistructstarts(n):bcistructends(n)) = tempbci(1:(bcistructends(n)-bcistructstarts(n)+1));
end
bcistruct.calibrationnev = calibrationnevstruct;
bcistruct.bcinev = bcinevstruct;
bcistruct.calibrationdata = trainingdat;
bcistruct.bcidata = bcidat;
bcistruct.modelparamsstruct = modelparamsstruct;
disp('saving struct; this could take a while...')
%save([filepath,bcifile{1},'_bcistruct'],'bcistruct')
end