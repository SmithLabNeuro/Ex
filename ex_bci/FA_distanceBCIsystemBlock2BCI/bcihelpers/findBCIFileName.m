function filename = findBCIFileName(filepath)
thismonkey = 'pe';%filepath(4:5); %folder has format D:monkeyname
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
                    userChoice = input([thisfilename,'.nev\nThis file? [y/n] ', ],'s');
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
            fprintf(['Using the following file: \n',filename,'.nev\n'])  
            fprintf('If you want a different file, enter the filename\n as the first input to this function (do not include .nev extention)\n')
        end
end