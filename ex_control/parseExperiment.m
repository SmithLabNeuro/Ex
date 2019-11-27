function m = parseExperiment(varargin) %#ok<*STOUT>

%Modified 03Feb2013 by ACS: fixed bug that caused punt when one of the
%input arguments was empty.

id = 1;
m = cell(0);
allSets = cell(0);

for s = 1:numel(varargin)
    
    allSets{s} = struct('pVal',{},'pName',{},'pDat',{},'pClass',{},'nDat',{});
    if isempty(varargin{s}), continue; end;
    
    c = varargin{s}.Children;
    pVal = cell(length(c),1);
    pName = cell(length(c),1);
    pDat = cell(length(c),1);
    pClass = cell(length(c),1);
    nDat = cell(length(c),1);
    
    for i = 1:length(c)
        if strcmp(c(i).Name,'grouped');
            for j = 1:length(c(i).Children)
                pName{i}{j} = c(i).Children(j).Name;
                pDat{i}{j} = eval(c(i).Children(j).Children.Data);
                pClass{i}{j} = class(pDat{i}{j}); %keep track of the class of the data to deal with strings properly -04Sep2013 ACS
                switch pClass{i}{j} %count the number of data for each parameter
                    case 'char'
                        pDat{i}{j} = pDat{i}{j}'; %transpose so that the 'words' are columns -04Sep2013 ACS
                        nDat{i}(j) = 1; %count a string as a single datum
                    otherwise
                        nDat{i}(j) = size(pDat{i}{j},2); %each column is a different datum, all items within a column go together (e.g., like params) -04Sep2013 ACS
                end;
            end
            if numel(unique(nDat{i}))>1, %adding this to avoid assuming that all grouped items were well-defined (have the same number of entries). -04Sep2013 ACS
                warning('readExperiment:unequalGrouping','Grouped items have different numbers of entries. Using the first %d entries for all items...',min(nDat{i}));
            end;
            pVal{i} = 1:min(nDat{i});
        else
            pName{i} = c(i).Name;
            pDat{i} = eval(c(i).Children(1).Data);
            pClass{i} = class(pDat{i}); %keep track of the class of the data to deal with strings properly -04Sep2013 ACS
            switch pClass{i} %count the number of data for each parameter
                case 'char'
                    pDat{i} = pDat{i}'; %transpose so that the 'words' are columns -04Sep2013 ACS
                    nDat{i} = 1; %count a string as a single datum
                otherwise
                    nDat{i} = size(pDat{i},2); %each column is a different datum, all items within a column go together (e.g., like params) -04Sep2013 ACS
            end;
            pVal{i} = 1:nDat{i};
        end
    end
    allSets{s} = struct('pVal',pVal,'pName',pName,'pDat',pDat,'pClass',pClass,'nDat',nDat);
end;

if numel(allSets)>1
    for s = numel(allSets):-1:2        
        pNamesUpdate = {allSets{s}.pName};
        if isempty(pNamesUpdate), continue, end; %If the update set is empty, then nothing needs done. -ACS 03Feb2014
        pNamesDefault = {allSets{s-1}.pName};
        groupedVarsDefault = cellfun(@iscell,pNamesUpdate);
        groupedVarsUpdate = cellfun(@iscell,pNamesUpdate);        
        assert(sum(groupedVarsDefault)<=1&&sum(groupedVarsUpdate<=1),'Use only one set of grouped variables!'); %this is too complicated to deal with right now if there is more than one set of grouped vars. -acs24jun2016
        %remove grouped variables from the list (cells don't work with
        %ismember):
        pNamesDefault(groupedVarsDefault) = {'groupedVar-doNotUpdate'};
        pNamesUpdate(groupedVarsUpdate) = {'groupedVar-doNotUse'};
        if any(groupedVarsDefault),
            %update default grouped vars
            if any(groupedVarsUpdate),
                allSets{s-1}(groupedVarsDefault) = allSets{s}(groupedVarsUpdate);
            end;
        elseif any(groupedVarsUpdate),
            %append the new group to the list...
            allSets{s-1} = vertcat(allSets{s-1},allSets{s}(groupedVarsUpdate));
        end;
        overwrite = find(ismember(pNamesDefault,pNamesUpdate)); %identify which values in the default set are overwritten by the higher-priority set
        for ow = 1:numel(overwrite)             
            allSets{s-1}(overwrite(ow)) = allSets{s}(ismember(pNamesUpdate,pNamesDefault(overwrite(ow)))); %replace the values to be overwritten with the new values
        end;     
        %now match the grouped names so they aren't appended...
        %-acs24jun2016 (sorry)
        pNamesDefault(groupedVarsDefault) = {'groupedVar'};
        pNamesUpdate(groupedVarsUpdate) = {'groupedVar'};
        allSets{s-1} = vertcat(allSets{s-1},allSets{s}(~ismember(pNamesUpdate,pNamesDefault))); %append new values to default set that are not present
    end;
end;

%put the structure back into the multiple cell array variables:
allSets = allSets{1};
fnames = fieldnames(allSets);
allSets = struct2cell(allSets); %#ok<NASGU>
for fn = 1:numel(fnames)
    eval(sprintf('%s=allSets(fn,:);',fnames{fn}));
end;

ndGridCmd = '[';
for i = 1:length(pVal)
    ndGridCmd = sprintf('%sp{%i} ',ndGridCmd,i);
end
ndGridCmd = [ndGridCmd(1:end-1) '] = ndgrid('];
for i = 1:length(pVal)
    ndGridCmd = sprintf('%spVal{%i},',ndGridCmd,i);
end
ndGridCmd = [ndGridCmd(1:end-1) ',1);'];

p = cell(length(pVal),1);
eval(ndGridCmd);

for i = 1:length(p{1}(:))
    
    for s = 1:numel(varargin) %go in reverse order of priority in case anything needs to get overwritten
        for sub = 1:numel(varargin{s})
            for j = 1:length(varargin{s}(sub).Attributes)
                eval(sprintf('m{id}.%s = varargin{s}(sub).Attributes(j).Value;',varargin{s}(sub).Attributes(j).Name));
            end
        end;
    end;
    
    for j = 1:length(pVal)
        if iscell(pName{j}) %i.e., grouped. Only really makes sense for conditions.
            for k = 1:length(pName{j})
                % WORK HERE
                eval(sprintf('m{id}.%s = pDat{j}{k}(:,p{j}(i))'';',pName{j}{k})); %added colon operator so columns are kept together (Each column is a separate condition), and transpose operator to turn columns into row vectors (and make words again for strings) -04Sep2013 ACS
            end
        else
            eval(sprintf('m{id}.%s = pDat{j}(:,p{j}(i))'';',pName{j})); %added colon operator so columns are kept together (Each column is a separate condition), and transpose operator to turn columns into row vectors (and make words again for strings) -04Sep2013 ACS
        end
    end
    id = id+1;
end