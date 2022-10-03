function updateExperimentSessionInDatabase(sqlDb, sessionNumber, subject, varargin)

if isempty(sqlDb)
    error('No database linked, so nothing written to database.')
end
if ~isempty(varargin)
    % write out the update string for the experiment_session with the given
    % sessionNumber
    updateStrCell = cell(1,length(varargin)/2);
    for fieldToUpdateInd = 1:2:length(varargin)
        fieldToUpdate = varargin{fieldToUpdateInd};
        valueToPutIn = varargin{fieldToUpdateInd + 1};
        % updating in SQL just involves setting the field to the value (and
        % it kind of ignores whether it's quoted or not, and relies on what
        % the DB says is the field type; quoting is a good generalization
        % to allow for both numbers and stringss
        % NOTE that this means your values all need to be strings, though
        % >.> (even if updating an int field, num2str the value before
        % running this function)
        updateStr = sprintf('%s = "%s"', fieldToUpdate, valueToPutIn);
        updateStrCell{(fieldToUpdateInd+1)/2} = updateStr;
    end
    % can make a bunch of updates by linking them with 'and's.
    updateStrTotal = strjoin(updateStrCell, ', ');
    
    % execute the database update
    sqlDb.exec(sprintf('UPDATE experiment_session SET %s WHERE session_number = %d AND animal="%s"', updateStrTotal, sessionNumber, subject));
end