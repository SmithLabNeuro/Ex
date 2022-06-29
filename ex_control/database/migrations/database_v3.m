% adding water_ml and threshold_rms fields to experiment_session table to
% track more session info in the database
exGlobals;
localDataDir = params.localDataDir;
homedirCont = dir('~');
homedir = homedirCont(1).folder;
homeTag = find(localDataDir == '~');
localDataDir = [localDataDir(1:homeTag-1) homedir localDataDir(homeTag+1:end)];
sqlDbPath = fullfile(localDataDir, 'database', 'experimentInfo.db');
sqlDb = sqlite(sqlDbPath);

%% update experiment_session to add water_ml and threshold_rms fields
% change move over experiment_session to experiment_session_old
alterStatement = 'ALTER TABLE experiment_session RENAME TO experiment_session_old';
sqlDb.exec(alterStatement);

% create experiment_session table
createStatement = [
    'CREATE TABLE experiment_session '...
        '(session_number INT, '...
        'date DATETIME, '...
        'animal CHAR, '...
        'rig CHAR, '...
        'notes TEXT NULLABLE, '...
        'water_ml INT, '... new field
        'threshold_rms INT, '... new field
        'FOREIGN KEY (rig) REFERENCES rig (name) '...
            'ON DELETE CASCADE, '...
        'FOREIGN KEY (animal) REFERENCES animal (name) '...
            'ON DELETE CASCADE, '...
        'CONSTRAINT session_pk PRIMARY KEY (session_number, animal))'
     ];
sqlDb.exec(createStatement);

% move data into this new table that includes water and threshold info for
% the session
insertStatement = ['INSERT INTO experiment_session '...
    '(session_number, date, animal, rig, notes)'...
    'SELECT session_number, date, animal, rig, notes FROM experiment_session_old'];
sqlDb.exec(insertStatement);

% drop old table
deleteStatement = 'DROP TABLE experiment_session_old;';
sqlDb.exec(deleteStatement);

