function [globals] = readUserXML_rig(stem)

% function [m params randoms globals] = readUserXML(subject,xmlFile)
%
% Basically a modified version of readExperiment that reads
% subject-specific settings. -ACS 05Sep2014

% reads an xml file and returns the 3 sets of parameters:
% conditions: parameters that are tied to condition number
% params: parameters that are fixed
% randoms: parameters that vary randomly on every trial, but are not tied
%   to condition

globals = struct([]);

userFile = ['rig_' lower(stem) '.xml'];

if ~exist(userFile,'file')
    return;
end

x = xml2struct(userFile);
x = removeComments(x);
x = removeBlanks(x);


for kid = 1:length(x.Children)
    node = x.Children(kid);
    if strcmp(node.Name,'globals')
        globals = node;
    else
        error('Not Valid:rig has field other than globals');
    end
end

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
