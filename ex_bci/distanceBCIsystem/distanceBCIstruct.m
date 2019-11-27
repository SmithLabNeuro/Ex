copyToHive
%filepath ='/Volumes/Queen/data/pepelepew/'; % for MAC
filepath = 'Z:data\pepelepew\'; % for windows
localPath = 'D:pepelepew\';
%filestems = {'Pe180928_s449a_distanceStabilityBCI'};
% 
% 
% 
filestems = [{'Pe181009_s457a_distanceStabilityBCI'}];


addpath('../distinguishable_colors/');
%error('ADD PATH TO NEURAL_NET_SORT_STUFF')
set(groot,'defaultAxesColorOrder',distinguishable_colors(16))

for fileinds = 1:length(filestems)  
    filestem = filestems{fileinds};
    calfilename =[filestem,'calib_0001'];
    bcifilename =[{[filestem,'_0002']}];
    hmmpostname = '';
%     if fileinds ~= 4
%         bcifilename = [filestem,'0002'];
%     else
%         bcifilename = [filestem,'0003'];
%     end
    parameterfile = [{[calfilename,'_lowdHmm0',hmmpostname]}];
    logfile = [{[calfilename,'_lowdHmm0',hmmpostname,'_log1.txt']}];
    addpath(genpath('../'))
    %load([filepath,bcifilename{1},'_bcistruct'])
    bcistruct = createBCIStruct(localPath,filepath,calfilename,bcifilename,logfile,parameterfile);
   % bcistruct.calibrationdata = getNS5Data(bcistruct.calibrationdata ,filepath,calfilename,'.ns5');
   % bcistruct.bcidata = getNS5Data(bcistruct.bcidata ,filepath,bcifilename,'.ns5');
    bcistruct.bcidata = removeBadTrials(bcistruct.bcidata);
    save([filepath,bcifilename{end},'_bcistruct.mat'],'bcistruct')
   % filedate = strsplit(filestem,'_');
   % newsavepath = ['C:/Users/smithlab/Dropbox/smithlabdata/bcidailyplots/',filedate{1}];
   % dailyDiagnostic(bcistruct,[bcifilename{end},'_bcistruct.mat'],newsavepath,true,2)
   % otherDailyPlots(bcistruct,[bcifilename{end},'_bcistruct.mat'],newsavepath)
  %  aggregatePlots('180813',filestem(3:8),newsavepath,'savePlots',true,'filePath',filepath);
end


