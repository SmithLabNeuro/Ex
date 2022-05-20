% changing up some primary key associations (to what they should have been
% to begin with >.>)
exGlobals;
localDataDir = params.localDataDir;
homedirCont = dir('~');
homedir = homedirCont(1).folder;
homeTag = find(localDataDir == '~');
localDataDir = [localDataDir(1:homeTag-1) homedir localDataDir(homeTag+1:end)];
sqlDbPath = fullfile(localDataDir, 'database', 'experimentInfo.db');
sqlDb = sqlite(sqlDbPath);

%% start with the primary key setup for experiment_session
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
        'FOREIGN KEY (rig) REFERENCES rig (name) '...
            'ON DELETE CASCADE, '...
        'FOREIGN KEY (animal) REFERENCES animal (name) '...
            'ON DELETE CASCADE, '...
        'CONSTRAINT session_pk PRIMARY KEY (session_number, animal))'
     ];
sqlDb.exec(createStatement);

% move data into this new, primary-keyed table
insertStatement = 'INSERT INTO experiment_session SELECT * FROM experiment_session_old';
sqlDb.exec(insertStatement);


%% now we do task_info primary key...
% change move over task_info to task_info_old
alterStatement = 'ALTER TABLE task_info RENAME TO task_info_old';
sqlDb.exec(alterStatement);

% create task_info table
createStatement = 'CREATE TABLE task_info ( name char PRIMARY KEY, parent_task char NULLABLE)';
sqlDb.exec(createStatement);

% move data into this new, primary-keyed table
insertStatement = 'INSERT INTO task_info SELECT * FROM task_info_old';
sqlDb.exec(insertStatement);

%% now we can do experiment_info because its relations have been updated
% change move over experiment_info to experiment_info_old
alterStatement = 'ALTER TABLE experiment_info RENAME TO experiment_info_old';
sqlDb.exec(alterStatement);

% create experiment_info table
createStatement = [
    'CREATE TABLE experiment_info '...
        '( session INT, '...
        'animal CHAR, '...
        'start_time DATETIME PRIMARY KEY, '...
        'task CHAR, '...
        'behavior_output_name CHAR, '...
        'neural_output_name CHAR, '...
        'experiment_results TEXT, '...
        'FOREIGN KEY (session, animal) REFERENCES experiment_session (session_number, animal) '...
            'ON DELETE CASCADE, '...
        'FOREIGN KEY (task) REFERENCES task_info (name) '...
        'ON DELETE CASCADE )'
     ];
sqlDb.exec(createStatement);

% move data into this new, primary-keyed table
insertStatement = ['INSERT INTO experiment_info '...
    '(session, start_time, task, behavior_output_name, neural_output_name, experiment_results) '...
    'SELECT session, start_time, task, behavior_output_name, neural_output_name, experiment_results FROM experiment_info_old'];
sqlDb.exec(insertStatement);
updateStatement = 'UPDATE experiment_info SET animal=''satchel''';
sqlDb.exec(updateStatement);


deleteStatement = 'DROP TABLE experiment_info_old;';
sqlDb.exec(deleteStatement);
deleteStatement = 'DROP TABLE experiment_session_old;';
sqlDb.exec(deleteStatement);
deleteStatement = 'DROP TABLE task_info_old;';
sqlDb.exec(deleteStatement);

