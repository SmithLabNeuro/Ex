function [decoded, times] = exDecode(codeIn)
%To decode ex codes into something resembling English. adam@adamcsnyder.com
%23Apr2012

%removed references to "ex_codes", which are now out of date. Anything that
%had been relying on "ex_codes" needs to be updated to use exGlobals
%-acs18may2016

%note also that there are two nearly identical versions of this function
%with the same name... but I'm averse to removing one because the rigs
%don't have access to utils, and people's workstations don't have access to
%smithlabrig... -acs18may2016
codes = [];
    exGlobals;   
    times = [];
    if isa(codeIn,'struct')&&isfield(codeIn,'codes') %added to make it easier to work with allCodes variables output by Ex -ACS 24AUG2012
        codeIn = codeIn.codes;
    end;
    if sum(size(codeIn)>1)>1
        if size(codeIn,2)==3, %indicates this is output from readNEV
            digCodes = codeIn(:,1)==0;
            exCodes = ismember(codeIn(:,2),1:255);
            codeIn = codeIn(digCodes&exCodes,2:3);
        end;
        times = codeIn(:,2);
    end;
    codeIn = codeIn(:,1);    
    codeNames = fieldnames(codes);
    codeNums = cell2mat(struct2cell(codes));
    nameInds = cell2mat(cellfun(@find,cellfun(@eq,repmat({codeNums},size(codeIn)),num2cell(codeIn),'uniformoutput',0),'uniformoutput',0));
    try
        times = times(ismember(codeIn,codeNums));
    catch %#ok<CTCH>
        times = [];
    end;
    decoded = codeNames(nameInds);
