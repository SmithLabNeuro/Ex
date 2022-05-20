function [sessionInfo, sessionNotes] = writeExperimentSessionToDatabase(sqlDb, params)

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
    'FROM experiment_info JOIN experiment_session '... we need some stuff from experiment_session as well (namely, the animal), so we're joining tables here
    'ON experiment_info.session=experiment_session.session_number '... joining on the session being correct
    'WHERE datetime(start_time) >= datetime(''now'', ''-1 day'', ''localtime'') AND '... only grabbing sessions in the past 24 hours
    'experiment_session.animal=''%s'''], params.SubjectID)); % grabbing sessions for the specific animal

% There is already some experiment info which might indicate a session has
% begun
if ~isempty(infoPossiblyRelated)
    % Looking at how many hours it has been since the task was run
    hoursFromStart = cell2mat(infoPossiblyRelated(:, 2));
    [sortedHours, sortInd] = sort(hoursFromStart);
    infoPossiblyRelated = infoPossiblyRelated(sortInd, :);
    if any(sortedHours<12)
        % Assume anything within 12 hours is the same session (is this
        % appropriate? I mean... unless you run a monkey twice in one day
        % and consider it a different session...)
        infoId = infoPossiblyRelated{sortedHours<12, 1};
        sessionAllInfo = sqlDb.fetch(sprintf('SELECT session_number, ifnull(notes,"") FROM experiment_session JOIN experiment_info ON experiment_info.session = experiment_session.session_number WHERE experiment_info.rowid = %d AND experiment_session.animal=''%s''', infoId, params.SubjectID));
        sessionInfo = sessionAllInfo{1};
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
        sqlDb.insert('experiment_session', {'session_number', 'date', 'animal', 'rig'}, {newSessionNumber, datestr(today, 'yyyy-mm-dd'), params.SubjectID, params.machine})
        sessionInfo = sqlDb.fetch(sprintf('SELECT session_number FROM experiment_session WHERE experiment_session.animal=''%s'' ORDER BY session_number DESC LIMIT 1', params.SubjectID));
        sessionInfo = sessionInfo{1};
        sessionNotes = [];
    end
else
    % no experiments within 1 day, but there still might be an
    % experiment_session *from today* which we check for here
    todaySessionAllInfo = sqlDb.fetch(sprintf('SELECT session_number, ifnull(notes,"") FROM experiment_session WHERE strftime(''%%Y-%%m-%%d'', date) == strftime(''%%Y-%%m-%%d'', ''now'', ''localtime'') and animal=''%s''', params.SubjectID));
    if size(todaySessionAllInfo, 1)==1
        % this might happen if runex was run but with no
        % experiment--here we avoid writing another row to
        % experiment_session
        sessionInfo = todaySessionAllInfo{1};
        sessionNotes = todaySessionAllInfo{2};
    elseif size(todaySessionAllInfo, 1)>1
        keyboard
    else
        sqlDb.insert('experiment_session', {'session_number', 'date', 'animal', 'rig'}, {newSessionNumber, datestr(today, 'yyyy-mm-dd'), params.SubjectID, params.machine})
        sessionInfo = sqlDb.fetch(sprintf('SELECT session_number FROM experiment_session WHERE animal=''%s'' ORDER BY session_number DESC LIMIT 1', params.SubjectID));
        sessionInfo = sessionInfo{1};
        sessionNotes = [];
    end
end