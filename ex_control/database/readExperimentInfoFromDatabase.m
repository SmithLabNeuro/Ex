function dbResult = readExperimentInfoFromDatabase(outfilename, varargin)
% reads experiment info entry
% - outfilename used to identify experiment_info being read
% varargin tells us what fields read, and if empty returns everything
global sqlDb 
if ~isempty(varargin)
    % because Matlab's (2019b) sqlite is dumb, we need to deal with
    % potentially null fields by explicitly returning empty strings...
    readFieldStr = [' ifnull(', strjoin(varargin, ', ""), ifnull('), ', "") '];
    dbResult = sqlDb.fetch(sprintf('SELECT %s FROM experiment_info WHERE behavior_output_name = "%s"', readFieldStr, outfilename));
else
    dbResult = sqlDb.fetch(sprintf('SELECT * FROM experiment_info WHERE behavior_output_name = "%s"', outfilename));
end
