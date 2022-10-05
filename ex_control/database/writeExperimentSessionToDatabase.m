function [sessionNumber, sessionNotes] = writeExperimentSessionToDatabase(sqlDb, params)

if isempty(sqlDb)
    error('No database linked, so nothing written to database.')
end

% grabbing the previous session number
lastSessionNumber = sqlDb.fetch(sprintf('SELECT ifnull(max(session_number), 0) from experiment_session where animal=''%s''', params.SubjectID));

% Setting a new session number, if needed
newSessionNumber = lastSessionNumber{1} + 1;

% Checking if there is any info that could have been from today... (this
% kind of takes into account someone running a monkey at 11PM and then a
% task at 1AM... so it might be a different calendar day but same actual
% day. Should this ever happen? No. Would this ever happen? When it comes
% to dates and times on computers who knows.)
infoPossiblyRelated = sqlDb.fetch(sprintf(['SELECT experiment_info.rowid, '... rowid from experiment_info (not experiment_session) to help identify row
    '(strftime(''%%s'', ''now'', ''localtime'') - strftime(''%%s'',start_time))/3600 '... find the number of HOURS between now and the task starts
    'FROM experiment_info LEFT OUTER JOIN experiment_session '... we need some stuff from experiment_session as well (namely, the animal), so we're joining tables here
    'ON experiment_info.session=experiment_session.session_number '... joining on the session being correct
    'WHERE datetime(start_time) >= datetime(''now'', ''-1 day'', ''localtime'') AND '... only grabbing sessions in the past 24 hours
    'experiment_info.animal=''%s'''], params.SubjectID)); % grabbing sessions for the specific animal

% here's to check if, even if there are no experiments within 1 day, there
% was an experiment_session *from today* (could happen if we run runex
% without recording any data; the session is still written out to allow for
% not taking, etc.)
todaySessionAllInfo = sqlDb.fetch(sprintf('SELECT session_number, ifnull(notes,"") FROM experiment_session WHERE strftime(''%%Y-%%m-%%d'', date) == strftime(''%%Y-%%m-%%d'', ''now'', ''localtime'') and animal=''%s''', params.SubjectID));


if size(todaySessionAllInfo, 1)==1
    % this might happen if runex was run but with no
    % experiment--here we avoid writing another row to
    % experiment_session
    sessionNumber = todaySessionAllInfo{1};
    sessionNotes = todaySessionAllInfo{2};
elseif size(todaySessionAllInfo, 1)>1
    % why would there be multiple sessions from today?
    keyboard
elseif ~isempty(infoPossiblyRelated)
    % no experiment_session entries with today's date, but there are
    % experiment_infos which might indicate a session from "today"--i.e.
    % within some hours but there was a day change
    
    % Looking at how many hours it has been since the task was run
    hoursFromStart = cell2mat(infoPossiblyRelated(:, 2));
    [sortedHours, sortInd] = sort(hoursFromStart);
    infoPossiblyRelated = infoPossiblyRelated(sortInd, :);
    if any(sortedHours<12)
        % Assume anything within 12 hours is the same session (is this
        % appropriate? I mean... unless you run a monkey twice in one day
        % and consider it a different session...)
        infoId = infoPossiblyRelated{sortedHours<12, 1};
        % 'SELECT session_number, ifnull(notes,"") FROM experiment_session JOIN experiment_info ON experiment_info.session = experiment_session.session_number WHERE experiment_info.rowid = 1958 AND experiment_session.animal='satchel''
        sessionAllInfo = sqlDb.fetch(sprintf('SELECT session_number, ifnull(notes,"") FROM experiment_session LEFT OUTER JOIN experiment_info ON experiment_info.session = experiment_session.session_number WHERE experiment_info.rowid = %d AND experiment_info.animal=''%s''', infoId, params.SubjectID));
        sessionNumber = sessionAllInfo{1};
        sessionNotes = sessionAllInfo{2};
    elseif any(sortedHours<16)
        % Not really sure why this exact spacing might happen--maybe if you
        % run a monkey at 6PM one day and start up at 9AM the next, but
        % oof.
        keyboard
    else
        % there was an experiment within 1 day (could be up to
        % 47.9999 hours, since the 'day' metric could be day 1 at
        % 12:01 AM vs day 2 at 11:59 PM, but that's still day 2 -
        % day 1 = 1 day...), but >16 hours so not counted as the same
        % session
        %
        % grab the rig name from the database to link to the experiment session
        rigName = sqlDb.fetch(sprintf('SELECT name FROM rig WHERE control_computer_name="%s"', params.machine));
        rigName = rigName{1};
        
        sqlDb.insert('experiment_session', {'session_number', 'date', 'animal', 'experimenter', 'rig'}, {newSessionNumber, datestr(today, 'yyyy-mm-dd'), params.SubjectID, params.experimenter, rigName})
        sessionNumber = sqlDb.fetch(sprintf('SELECT session_number FROM experiment_session WHERE experiment_session.animal=''%s'' ORDER BY session_number DESC LIMIT 1', params.SubjectID));
        sessionNumber = sessionNumber{1};
        sessionNotes = [];
    end
else
    % grab the rig name from the database to link to the experiment session
    rigName = sqlDb.fetch(sprintf('SELECT name FROM rig WHERE control_computer_name="%s"', params.machine));
    rigName = rigName{1};
    sqlDb.insert('experiment_session', {'session_number', 'date', 'animal', 'experimenter', 'rig'}, {newSessionNumber, datestr(today, 'yyyy-mm-dd'), params.SubjectID, params.experimenter, rigName})
    sessionNumber = sqlDb.fetch(sprintf('SELECT session_number FROM experiment_session WHERE animal=''%s'' ORDER BY session_number DESC LIMIT 1', params.SubjectID));
    sessionNumber = sessionNumber{1};
    sessionNotes = [];
end