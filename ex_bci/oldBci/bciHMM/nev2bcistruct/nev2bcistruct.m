function [ dat ] = nev2bcistruct(filename,varargin)
% nev_info is an optional argument

addpath helpers
% optional args
nevreadflag = 0;
nev_info = [];
assignopts (who, varargin)
%% important codes
starttrial = 1;
endtrial = 255;

if nevreadflag ==1
    nev = filename;
else
    if isempty(findstr(filename,'.nev'))
        nev = readNEV([filename,'.nev']);
        nev_info = NEV_displayheader([filename,'.nev']);
    else
        nev = readNEV(filename);
        nev_info = NEV_displayheader(filename);
    end
end
digcodes = nev(nev(:,1)==0,:);
%spikecodes = nev(nev(:,1)~=0,:);
diginnevind = find(nev(:,1)==0);
trialstartindstemp = (find(digcodes(:,2)==starttrial));
trialstartinds = diginnevind(trialstartindstemp);
trialstarts = nev(trialstartinds,3);

trialendindstemp = (find(digcodes(:,2)==endtrial));
trialendinds = diginnevind(trialendindstemp);
trialends = nev(trialendinds,3);
[trialstarts, trialends,trialstartgood,trialendgood] = detectMissingStartEndCode(trialstarts,trialends);
trialstartinds = trialstartinds(trialstartgood);
trialendinds = trialendinds(trialendgood);

if length(trialstarts)~=length(trialends) || sum((trialends-trialstarts)<0)
    % fix it
    if sum(trialstarts(1:end-1)>=trialends)==0
        trialstarts = trialstarts(1:end-1);
    end
end

%% get session initial params
predatcodes = digcodes(digcodes(:,3)<trialstarts(1),:);
tempdata.text = char(predatcodes(predatcodes(:,2)>=256 & predatcodes(:,2)<512,2)-256)';
tempdata = getDatParams(tempdata);

%% Make Struct
%spikecodetimes = spikecodes(:,3);
%digcodetimes = digcodes(:,3);
for n = 1:length(trialstarts)
    if mod(n,100) == 0
        fprintf('Processed nev for %i trials of %i...\n',n,length(trialstarts))
    end
    dat(n).time = [trialstarts(n) trialends(n)];
    thisnev = nev(trialstartinds(n):trialendinds(n),:);
    trialdig = thisnev(thisnev(:,1)==0,:);
    tempspikes = thisnev(thisnev(:,1)~=0,:);
    tempspikes(:,3) = tempspikes(:,3)*30000;
    %digind1 = binaryTimeSearch(digcodetimes,trialstarts(n));
    %digind2 = binaryTimeSearch(digcodetimes,trialends(n));
    
   % indices = binaryTimeSearch(digcodetimes,trialstarts(n),trialends(n));
    %trialdig = digcodes(digcodetimes>=trialstarts(n) & digcodetimes<=trialends(n),:);
    %trialdig = digcodes(digind1:digind2,:);
    dat(n).text = char(trialdig(trialdig(:,2)>=256 & trialdig(:,2)<512,2)-256)';

    dat(n).trialcodes = trialdig(trialdig(:,2)<256,:);
    trialdig(:,3) = trialdig(:,3)*30000;
    dat(n).events = uint32(trialdig);
    %indices = binaryTimeSearch(spikecodetimes,trialstarts(n),trialends(n));
    % indices = spikecodetimes>=trialstarts(n) & spikecodetimes<=trialends(n);
    %tempspikes = spikecodes(indices,:);
    %tempspikes(:,3) = tempspikes(:,3)*30000;
    dat(n).firstspike = tempspikes(1,3);
    dat(n).spiketimesdiff = uint16(diff(tempspikes(:,3)));
    dat(n).spikeinfo = uint16(tempspikes(:,1:2));
    dat(n).result = dat(n).events(dat(n).events(:,2)>=150 & dat(n).events(:,2)<=158,2);
    dat(n).params.block = tempdata.params.trial;
    if ~isempty(nev_info); dat(n).nevinfo.nevclockstart = nev_info.nevclockstart; end
    if(isempty(dat(n).result))
        dat(n).result = NaN;
    end
end
dat = getDatParams(dat);

end