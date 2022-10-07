% correctly referencing new experiment_session from experiment_info
%   this happened because when primary keys were added to
%   experiment_session, the first step was to move experiment_session to
%   experiment_session_old, but all the foreign keys from experiment_info
%   cascaded into the change--so now we need to relink to the new, correct
%   experiment_session
exGlobals;
localDataDir = params.localDataDir;
homedirCont = dir('~');
homedir = homedirCont(1).folder;
homeTag = find(localDataDir == '~');
localDataDir = [localDataDir(1:homeTag-1) homedir localDataDir(homeTag+1:end)];
sqlDbPath = fullfile(localDataDir, 'database', 'experimentInfo.db');
% we're gonna backup the database to the millisecond, to super minimize
% chances of someone accidentally double running this code and losing data
sqlDbBackupPath = fullfile(localDataDir, 'database', sprintf('experimentInfo_%s.db', datestr(now, 'yyyymmdd_HHMMSSFFF')));
sqlDb = sqlite(sqlDbPath);

%% add field to experiment_info
% change move over experiment_info to experiment_info_old
alterStatement = 'ALTER TABLE experiment_info RENAME TO experiment_info_old';
sqlDb.exec(alterStatement);

% create experiment_info table with new field
createStatement = [
    'CREATE TABLE experiment_info '...
        '( session INT, '...
        'animal CHAR, '...
        'start_time DATETIME PRIMARY KEY, '...
        'task CHAR, '...
        'behavior_output_name CHAR, '...
        'neural_output_name CHAR, '...
        'experiment_results TEXT, '...
        'extra_notes TEXT, '... new field being added
        'FOREIGN KEY (session, animal) REFERENCES experiment_session (session_number, animal) '...
            'ON DELETE CASCADE, '...
        'FOREIGN KEY (task) REFERENCES task_info (name) '...
        'ON DELETE CASCADE )'
     ];
sqlDb.exec(createStatement);

% move data into this new, primary-keyed table
insertStatement = ['INSERT INTO experiment_info '...
    '(session, animal, start_time, task, behavior_output_name, neural_output_name, experiment_results, extra_notes) '...
    'SELECT session, animal, start_time, task, behavior_output_name, neural_output_name, experiment_results, extra_notes FROM experiment_info_old'];
sqlDb.exec(insertStatement);

deleteStatement = 'DROP TABLE experiment_info_old;';
sqlDb.exec(deleteStatement);

