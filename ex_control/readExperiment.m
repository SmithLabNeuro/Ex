function [cnds, params, randoms, globals] = readExperiment(xmlFile,subject,machine)
% function [m params randoms] = readExperiment(xmlFile)
%
% reads an xml file and returns the 3 sets of parameters:
% cnds: main attribute
% params: attribute, conditions and parameters
% randoms: parameters that vary randomly on every trial

cnds = struct();

x = xml2struct(xmlFile);
x = removeComments(x);
x = removeBlanks(x);


% name,repeats,ex,bgcolor must in stimulus
mainAtt = squeeze(struct2cell(x.Attributes)); %added squeeze -04Sep2013 ACS
if (isempty(find(strcmp(mainAtt,'name'), 1)) || isempty(find(strcmp(mainAtt,'repeats'), 1)) || ...
        isempty(find(strcmp(mainAtt,'ex'), 1)) || isempty(find(strcmp(mainAtt,'bgColor'), 1)))
    error('Missing main attributes');
end

params = cell(0);

% read in subject file
subject_params = struct([]);
subject_randoms = struct([]);
subject_globals = struct([]);
if nargin>1 %added to support subject-specific settings -05Sep2013 ACS
    [subject_params, subject_randoms, subject_globals] = readUserXML(subject,xmlFile);
    if ~isempty(fieldnames(subject_globals))
        subject_globals_value = squeeze(struct2cell(subject_globals));
        subject_globals_name = fieldnames(subject_globals);
        for i = 1:length(subject_globals_name)
            if isnumeric(subject_globals_value{i})
                if size(subject_globals_value{i},2)>1
                    error('Not Valid:global %s in subject should not changeable',subject_globals_name);
                else
                    subject_globals.(subject_globals_name{i}) = subject_globals_value{i}';
                end
            end
        end
    end
end

rig_globals = struct([]);
% read in rig file
if nargin>2
    [rig_globals] = readUserXML_rig(machine);
    if isempty(fieldnames(rig_globals)) || isempty(rig_globals.Children) 
        rig_globals = [];
    else
        for i = 1:length(rig_globals.Children)
            if isnumeric(eval(rig_globals.Children(i).Children.Data))
                if size(eval(rig_globals.Children(i).Children.Data),2)>1
                    error('Not Valid:global %s in rig should not changeable',rig_globals.Children(i).Children.Name);
                end
            end
        end
        rig_globals=parseExperiment(rig_globals);
        rig_globals = rig_globals{1};
    end
end

% combine globals
if nargin > 2
    if ~isempty(rig_globals)
        rig_globals_name = fieldnames(rig_globals);
        for i = 1:length(rig_globals_name)
            subject_globals.(rig_globals_name{i})=rig_globals.(rig_globals_name{i});
        end
    end
end
globals = subject_globals;



% check for duplicate params/randoms names in stimulus file
allName={};
for kid = 1:length(x.Children)
    if strcmp(x.Children(kid).Name,'params')
        for i = 1:length(x.Children(kid).Children)
            if strcmp(x.Children(kid).Children(i).Name,'grouped')
                for j = 1:length(x.Children(kid).Children(i).Children)
                    allName = cat(2,allName,{x.Children(kid).Children(i).Children(j).Name});
                end
            else
                allName = cat(2,allName,{x.Children(kid).Children(i).Name});
            end
        end
    elseif strcmp(x.Children(kid).Name,'randoms')
        if ~isempty(x.Children(kid).Children)
            for i = 1:length(x.Children(kid).Children)
                randoms_value = eval(x.Children(kid).Children(i).Children.Data);
                if size(randoms_value,1) ~=1
                    error('randoms is not valid');
                end
            end
        end
        allName = cat(2,allName,{x.Children(kid).Children.Name});
    end
end

for i = 1:length(allName)
    for j = i+1:length(allName)
        if strcmp(allName{i},allName{j})
            error('Stimulus file specified the value of %s more than once',allName{i});
        end
    end
end

% check for duplicate params/randoms names in subject file
subjectName = {};
if ~isempty(subject_params) && isfield(subject_params,'Children') && ~isempty(subject_params.Children)
    for i = 1:length(subject_params.Children)
        if strcmp(subject_params.Children(i).Name,'grouped')
            for j = 1:length(subject_params.Children(i).Children)
                subjectName = cat(2,subjectName,{subject_params.Children(i).Children(j).Name});
            end
        else
            subjectName = cat(2,subjectName,{subject_params.Children(i).Name});
        end
    end
end

if ~isempty(fieldnames(subject_randoms))
    subject_randoms_value = squeeze(struct2cell(subject_randoms));
    for i = 1:length(subject_randoms_value)
        if size(subject_randoms_value{i},1) ~= 1
            error('subject random is not valid');
        end
    end
    
    subjectName = cat(2,subjectName,fieldnames(subject_randoms)');
end
if length(subjectName)>1
    for i = 1:length(subjectName)
        for j = i+1:length(subjectName)
            if strcmp(subjectName{i},subjectName{j})
                error('Subject file specified the value of %s more than once',subjectName{i});
            end
        end
    end
end

% globals not in params
if ~isempty(globals)
    globals_name = fieldnames(globals);
    for i = 1:length(globals_name)
        if ismember(globals_name{i},allName) || ismember(globals_name{i},subjectName)
            error('Globals %s can not overwrite existing params',globals_name{i});
        end
        %params.(globals_name{i}) = globals.(globals_name{i});
    end
end

paramIndex = find(strcmp({x.Children.Name},'params'));
randomIndex = find(strcmp({x.Children.Name},'randoms'));
userFile = ['subject_',lower(subject),'.xml'];

% name,repeats,ex can not be anywhere else
if ~isempty(x.Children(paramIndex).Attributes)
    paramsAttr = squeeze(struct2cell(x.Children(paramIndex).Attributes));
    if ismember('name',paramsAttr(1,:)') || ismember('repeats',paramsAttr(1,:)') || ismember('ex',paramsAttr(1,:)')
        error('name,repeats,ex can not be anywhere else');
    end
end
if ismember('name',allName) || ismember('repeats',allName) || ismember('ex',allName) || ...
        ismember('name',subjectName) ||  ismember('repeats',subjectName) ||  ismember('ex',subjectName) || ...
        ismember('name',globals_name) || ismember('repeats',globals_name) || ismember('ex',globals_name)
    error('name,repeats,ex can not be anywhere else');
end

% type,runline can not be anywhere else
if ismember('type',allName) || ismember('runline',allName) || ...
        ismember('type',subjectName) || ismember('runline',subjectName) || ...
        ismember('type',globals_name) || ismember('runline',globals_name)
    error('type,runline can not be anywhere else');
end

% currentBlock, currentCnd can not exist anywhere
if ismember('currentBlock',allName) || ismember('currentCnd',allName) || ...
        ismember('currentBlock',subjectName) || ismember('currentCnd',subjectName) || ...
        ismember('currentBlock',globals_name) || ismember('currentCnd',globals_name)
    error('currentBlock,currentCnd can not exist anywhere');
end

%params overwrite  %check group params size
if exist(userFile,'file')
    isGrouped = strcmp({x.Children(paramIndex).Children.Name},'grouped');
    if sum(isGrouped)>0 % group exist in basic file
        for idx = find(isGrouped)
            if length(x.Children(paramIndex).Children(idx).Children)>1
                for i = 1:length(x.Children(paramIndex).Children(idx).Children)
                    for j = i+1:length(x.Children(paramIndex).Children(idx).Children)
                        if size(eval(x.Children(paramIndex).Children(idx).Children(i).Children.Data),2)~=size(eval(x.Children(paramIndex).Children(idx).Children(j).Children.Data),2)
                            error('Grouped items %s in original xml file have different numbers of entries.',x.Children(paramIndex).Children(idx).Children(i).Children.Name);
                        end
                    end
                end
            end
        end
        if ~isempty(subject_params) || ~isempty(fieldnames(subject_randoms))
            error('Original xml file are not allowed to be overwrited when grouped exists.');
        end
    else %group not exist in basic file
        if ~isempty(subject_params) && isfield(subject_params,'Children') && ~isempty(subject_params.Children)
            for i = 1:length(subject_params.Children)
                if strcmp(subject_params.Children(i).Name,'grouped')
                    if length(subject_params.Children(i).Children)>1
                        for ii = 1:length(subject_params.Children(i).Children)
                            for jj = ii+1:length(subject_params.Children(i).Children)
                                if size(eval(subject_params.Children(i).Children(ii).Children.Data),2)~=size(eval(subject_params.Children(i).Children(jj).Children.Data),2)
                                    error('Grouped items %s in subject have different numbers of entries.',subject_params.Children(i).Children(ii).Children.Name);
                                end
                            end
                        end
                    end
                    for j = 1:length(subject_params.Children(i).Children)
                        paramNameIndex = find(strcmp({x.Children(paramIndex).Children.Name},subject_params.Children(i).Children(j).Name));
                        if ~isempty(randomIndex)
                            randomNameIndex = find(strcmp({x.Children(randomIndex).Children.Name},subject_params.Children(i).Children(j).Name));
                        else
                            randomNameIndex = [];
                        end
                        if ~isempty(paramNameIndex)
                            x.Children(paramIndex).Children(paramNameIndex)=[];
                        elseif ~isempty(randomNameIndex)
                            x.Children(randomIndex).Children(randomNameIndex)=[];
                        end
                    end
                else
                    paramNameIndex = find(strcmp({x.Children(paramIndex).Children.Name},subject_params.Children(i).Name));
                    if ~isempty(randomIndex)
                        randomNameIndex = find(strcmp({x.Children(randomIndex).Children.Name},subject_params.Children(i).Name));
                    else
                        randomNameIndex = [];
                    end
                    if ~isempty(paramNameIndex)
                        x.Children(paramIndex).Children(paramNameIndex)=[];
                    elseif ~isempty(randomNameIndex)
                        x.Children(randomIndex).Children(randomNameIndex)=[];
                    end
                end
            end
            x.Children(paramIndex).Children = [x.Children(paramIndex).Children subject_params.Children];
        end
        % randoms overwrite
        if ~isempty(subject_randoms)
            randoms_field = fieldnames(subject_randoms);
            for i = 1:length(randoms_field)
                paramNameIndex = find(strcmp({x.Children(paramIndex).Children.Name},randoms_field{i}));
                if ~isempty(randomIndex)
                    randomNameIndex = find(strcmp({x.Children(randomIndex).Children.Name},randoms_field{i}));
                else
                    randomNameIndex = [];
                end
                if ~isempty(paramNameIndex)
                    x.Children(paramIndex).Children(paramNameIndex)=[];
                elseif ~isempty(randomNameIndex)
                    x.Children(randomIndex).Children(randomNameIndex)=[];
                end
            end
        end
    end
end

temp_params = x.Children(paramIndex);
temp_conditions.Name = 'conditions';
temp_conditions.Attributes = [];
temp_conditions.Data = '';
temp_conditions.Children = [];
for i = length(temp_params.Children):-1:1
    if (strcmp(temp_params.Children(i).Name,'grouped') || size(str2num(temp_params.Children(i).Children.Data),2)>1)
        temp_conditions.Children = [temp_params.Children(i) temp_conditions.Children];
        temp_params.Children(i)=[];
    end
end
if ~isempty(temp_conditions.Children)
    cnds = parseExperiment(temp_conditions);
    params.emptyCnd = 0;
else
    cnds.emptyCnd = 1;
    cnds = {cnds};
end
if ~isempty(temp_params.Children)
    temp_params = parseExperiment(temp_params);
    temp_params = temp_params{1};
end
temp_params_name = fieldnames(temp_params);
if ~isempty(temp_params_name)
    for i = 1:length(temp_params_name)
        params.(temp_params_name{i}) = temp_params.(temp_params_name{i});
    end
end


params.name = cell2mat(mainAtt(2,strcmp('name',mainAtt(1,:))));
params.rpts = str2double(cell2mat(mainAtt(2,strcmp('repeats',mainAtt(1,:)))));
params.exFileName = cell2mat(mainAtt(2,strcmp('ex',mainAtt(1,:))));
if (isempty(find(strcmp(allName,'bgColor'), 1)) && isempty(find(strcmp(subjectName,'bgColor'), 1)))
    params.bgColor = cell2mat(mainAtt(2,strcmp('bgColor',mainAtt(1,:))));
end

if ~isempty(randomIndex)
    for i = 1:length(x.Children(randomIndex).Children)
        thisParam = x.Children(randomIndex).Children(i);
        paramName = thisParam.Name;
        paramData = thisParam.Children.Data;
        eval(sprintf('subject_randoms.%s = %s;',paramName,paramData));
    end
end
randoms = subject_randoms;

end

function x = removeBlanks(x)
keep = ones(length(x.Children),1);
for i = 1:length(x.Children)
    if strcmp(x.Children(i).Name,'#text')
        if all(isspace(x.Children(i).Data))
            keep(i) = 0;
        else
            error('xml file does not have proper format (outside)');
        end
    end
end
x.Children = x.Children(find(keep));
for i = 1:length(x.Children)
    keep = ones(length(x.Children(i).Children),1);
    for j = 1:length(x.Children(i).Children)
        if strcmp(x.Children(i).Children(j).Name,'#text')
            if all(isspace(x.Children(i).Children(j).Data))
                keep(j) = 0;
            else
                error('xml file does not have proper format(inside)');
            end
        end
    end
    x.Children(i).Children = x.Children(i).Children(find(keep));
end
end

function x = removeComments(x)
keep = ones(length(x.Children),1);
for i = 1:length(x.Children)
    if strcmp(x.Children(i).Name,'#comment')
        keep(i) = 0;
    else
        x.Children(i) = removeComments(x.Children(i));
    end
end

x.Children = x.Children(find(keep));
end