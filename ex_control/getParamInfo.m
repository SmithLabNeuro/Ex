function [pVal pName pDat pClass nDat] = getParamInfo(c)

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