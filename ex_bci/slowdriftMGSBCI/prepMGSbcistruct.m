    [~,filename,ext]=fileparts(allfiles(fileind).name);
    if strcmp(ext,'.nev')
        dat = nev2sortedStruct([filepath,filename,'.nev'],0.2,0);
        nevinfo = NEV_displayheader([filepath,filename,'.nev']);
        for n = 1:length(dat)
            dat(n).nevinfo.nevclockstart = nevinfo.nevclockstart;
        end
        dat = getNS5Data(dat,filepath,filename,'.ns5',3,3,3);       
        save([filepath,filename,'_dat'],'dat');
        
    end