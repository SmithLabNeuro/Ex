function [params,randoms,globals] = readUserXML(subject,xmlFile)

% function [m params randoms globals] = readUserXML(subject,xmlFile)
%
% Basically a modified version of readExperiment that reads
% subject-specific settings. -ACS 05Sep2014

% reads an xml file and returns the 3 sets of parameters:
% conditions: parameters that are tied to condition number
% params: parameters that are fixed
% randoms: parameters that vary randomly on every trial, but are not tied
%   to condition

randoms = struct();
globals = struct();
params = struct([]);

if isempty(subject)
    return;
end

userFile = ['subject_' lower(subject) '.xml'];

x = xml2struct(userFile);
x = removeComments(x);
x = removeBlanks(x);


for kid = 1:length(x.Children)
    node = x.Children(kid);
    if strcmp(node.Name,'globals');
        for i = 1:length(node.Children)
            thisParam = node.Children(i);
            paramName = thisParam.Name;
            paramData = thisParam.Children.Data;
            eval(sprintf('%s.%s = %s;',node.Name,paramName,paramData));
        end;
    elseif strcmp(node.Name,'experiment')
        attr = squeeze(struct2cell(node.Attributes)); %added squeeze -04Sep2013 ACS
        exptName = cell2mat(attr(2,strcmp('name',attr(1,:))));
        if ~strcmp(exptName,xmlFile), continue; end %only read the section that includes the information for the current experiment
        for i  = 1:length(node.Children)
            if strcmp(node.Children(i).Name,'params')
                params = node.Children(i);
            elseif strcmp(node.Children(i).Name,'randoms')
                thisParam = node.Children(i).Children;
                paramName = thisParam.Name;
                paramData = thisParam.Children.Data;
                eval(sprintf('%s.%s = %s;',node.Children(i).Name,paramName,paramData));
            end
        end
    end;
end;

end

function x = removeBlanks(x)
keep = ones(length(x.Children),1);
for i = 1:length(x.Children)
    if strcmp(x.Children(i).Name,'#text')
        if all(isspace(x.Children(i).Data))
            keep(i) = 0;
        else
            error('xml file does not have proper format (outside experiment)');
        end
    end
end
x.Children = x.Children(find(keep));
for i = 1:length(x.Children)
    if strcmp(x.Children(i).Name,'globals')
        keep = ones(length(x.Children(i).Children),1);
        for j = 1:length(x.Children(i).Children)
            if strcmp(x.Children(i).Children(j).Name,'#text')
                if all(isspace(x.Children(i).Children(j).Data))
                    keep(j) = 0;
                else
                    error('xml file does not have proper format(globals)');
                end
            end
        end
        x.Children(i).Children = x.Children(i).Children(find(keep));
    else
        keep = ones(length(x.Children(i).Children),1);
        for j = 1:length(x.Children(i).Children)
            if strcmp(x.Children(i).Children(j).Name,'#text')
                if all(isspace(x.Children(i).Children(j).Data))
                    keep(j) = 0;
                else
                    error(sprintf('xml file does not have proper format(%s)',x.Children(i).Attributes.Value));
                end
            end
        end
        x.Children(i).Children = x.Children(i).Children(find(keep));
        for j = 1:length(x.Children(i).Children)
            keep1 = ones(length(x.Children(i).Children(j).Children),1);
            for t = 1:length(x.Children(i).Children(j).Children)
                if strcmp(x.Children(i).Children(j).Children(t).Name,'#text')
                    if all(isspace(x.Children(i).Children(j).Children(t).Data))
                        keep1(t) = 0;
                    else
                        error(sprintf('xml file does not have proper format(%s)',x.Children(i).Attributes.Value));
                    end
                end
            end
            x.Children(i).Children(j).Children = x.Children(i).Children(j).Children(find(keep1));
        end
    end
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
