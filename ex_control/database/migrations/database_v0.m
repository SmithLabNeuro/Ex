% an original set of migration files poofed to nonexistance, so here is the
% setup of the current database, and then we'll migrate from it
exGlobals;
localDataDir = params.localDataDir;
homedirCont = dir('~');
homedir = homedirCont(1).folder;
homeTag = find(localDataDir == '~');
localDataDir = [localDataDir(1:homeTag-1) homedir localDataDir(homeTag+1:end)];
%sqlDbPath = fullfile(localDataDir, 'database', 'experimentInfo.db');
%destroyDbMaybe = questdlg('Doing this might destroy an existing database!! Continue?');
%switch destroyDbMaybe
%    case 'Yes'
%        sqlDb = sqlite(sqlDbPath, 'create');
%    otherwise
%        error('Not destroying anything here...')
%end

databaseDir = fullfile(localDataDir, 'database');
sqlDbPath = fullfile(databaseDir, 'experimentInfo.db');

% create db directory if its not there
if ~exist(databaseDir, 'dir')
    mkdir(databaseDir);
end

% create DB (or destroy old one and create new one)
if exist(sqlDbPath, 'file')
    destroyDbMaybe = questdlg(sprintf('A database file already exists at %s, so running this will destroy an existing database!! Continue?', sqlDbPath));
    switch destroyDbMaybe
        case 'Yes'
            sqlDb = sqlite(sqlDbPath, 'create');
        otherwise
            error('Not destroying anything here...')
    end
else
    sqlDb = sqlite(sqlDbPath, 'create');
end

% create animal table
createStatement = 'CREATE TABLE animal (name char PRIMARY KEY, birth_year INT)';
sqlDb.exec(createStatement)

% create brain_area table
createStatement = 'CREATE TABLE brain_area (name char PRIMARY KEY,long_name char )';
sqlDb.exec(createStatement)

% create experiment_info table
createStatement = [
    'CREATE TABLE experiment_info '...
        '( session INT, '...
        'start_time DATETIME, '...
        'task CHAR, '...
        'behavior_output_name CHAR, '...
        'neural_output_name CHAR, '...
        'experiment_results TEXT, '...
        'FOREIGN KEY (session) REFERENCES experiment_session (session_number) '...
            'ON DELETE CASCADE, '...
        'FOREIGN KEY (task) REFERENCES task_info (name) '...
        'ON DELETE CASCADE )'
     ];
sqlDb.exec(createStatement);

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
            'ON DELETE CASCADE )'
     ];
sqlDb.exec(createStatement);

% create implant_info table
createStatement = [
    'CREATE TABLE implant_info '...
        '( animal char, '...
        'brain_area char, '...
        'implant_date datenum, '...
        'implant_type text '...
            'CHECK( implant_type IN ("chamber", "array") ), '...
        'detail text, '...
        'FOREIGN KEY (animal) REFERENCES animal (name) '...
            'ON DELETE CASCADE, '...
        'FOREIGN KEY (brain_area) REFERENCES brain_area (name) '...
            'ON DELETE CASCADE )'
        ];
sqlDb.exec(createStatement);

% create rig table
createStatement = 'CREATE TABLE rig (name char PRIMARY KEY, room char, building char )';
sqlDb.exec(createStatement);

% create task_info table
createStatement = 'CREATE TABLE task_info ( name char, parent_task char NULLABLE)';
sqlDb.exec(createStatement);