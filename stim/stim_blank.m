function stim_blank(optstr,w,objID,arg)
%function stim_blank(optstr,w,objID,arg)
%
% showex helper function for 'blank' stim class
%
% Each helper function has to have the ability to do 3 things:
% (1) parse the input arguments from the 'set' command and precompute
% anything that is necessary for that stimulus
% (2) issue the display commands for that object
% (3) clean up after that object is displayed (not always necessary)

global objects;
global sv;

if strcmp(optstr,'setup')
    a = sscanf(arg,'%i');
    % arguments: (1) frameCount
    stimname = mfilename;
    objects{objID} = struct('type',stimname(6:end),'frame',0,'fc',a(1));
elseif strcmp(optstr,'display')
    % do nothing for a blank
elseif strcmp(optstr,'cleanup')
    % nothing necessary for this stim class
else
    error('Invalid option string passed into stim_*.m function');
end

