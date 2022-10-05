% add control_computer_name to link machine name to rig name
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
copyfile(sqlDbPath, sqlDbBackupPath)
sqlDb = sqlite(sqlDbPath);


%% update experiment_session to add water_ml and threshold_rms fields
% change move over experiment_session to experiment_session_old
alterStatement = 'ALTER TABLE rig ADD COLUMN control_computer_name';
sqlDb.exec(alterStatement);
