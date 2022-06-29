function dbResult = readExperimentInfoFromDatabase(outfilename, varargin)
% reads experiment info entry
% - outfilename used to identify experiment_info being read, 
%   - if it's a char, it's assumed to be the primary key--the
%     behavior_output_name
%   - if it's a double, it's assumed to be the rowid, which SQLite keeps as
%     an autoincremented value for each row
% varargin tells us what fields read, and if empty returns everything
global sqlDb 

if ischar(outfilename)
    whereClause = sprintf('WHERE behavior_output_name = "%s"', outfilename);
elseif isnumeric(outfilename)
    whereClause = sprintf('WHERE rowid = %d', outfilename);
end

if ~isempty(varargin)
    % because Matlab's (2019b) sqlite is dumb, we need to deal with
    % potentially null fields by explicitly returning empty strings...
    readFieldStr = [' ifnull(', strjoin(varargin, ', ""), ifnull('), ', "") '];
    dbResult = sqlDb.fetch(sprintf('SELECT %s FROM experiment_info %s', readFieldStr, whereClause));
else
    dbResult = sqlDb.fetch(sprintf('SELECT * FROM experiment_info %s', whereClause));
end
