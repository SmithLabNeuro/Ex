function writeExperimentInfoToDatabase(sessionInfo, xmlParams, outfilename, varargin)
% writes or updates experiment entries
% - sessionInfo and xmlParams only needed if we're inserting a new value,
% which is indicated by 'varargin' being empty
% - outfilename lets us locate existing experiment_info's
global sqlDb params
if ~isempty(varargin)
    % at this point, experiment_info matching these parameters should
    % already have been written, so we link it to the neural data being
    % recorded
    updateStrCell = {};
    for fieldToUpdateInd = 1:2:length(varargin)
        fieldToUpdate = varargin{fieldToUpdateInd};
        valueToPutIn = varargin{fieldToUpdateInd + 1};
        % this bit allows us to append to certain fields
        if contains(fieldToUpdate, {'neural_output_name', 'extra_notes'})
            fieldCurrent = sqlDb.fetch(sprintf('SELECT ifnull(%s, "") FROM experiment_info WHERE behavior_output_name="%s"', fieldToUpdate, outfilename));
            if ~isempty(fieldCurrent{1})
                valueToPutIn = sprintf("%s\n%s", fieldCurrent{1}, valueToPutIn);
            end
        end
        updateStr = sprintf('%s = "%s"', fieldToUpdate, valueToPutIn);
        updateStrCell{end+1} = updateStr;
    end
    updateStrTotal = strjoin(updateStrCell, ', ');
    sqlDb.exec(sprintf('UPDATE experiment_info SET %s WHERE behavior_output_name = "%s"', updateStrTotal, outfilename));
    %             sqlDb.insert('experiment_info', {'start_time', 'task', 'session', 'behavior_output_name', 'neural_output_name'}, {datestr(now, 'yyyy-mm-dd HH:MM:SS'), xmlBase, sessionInfo, outfilename, neuralOutName})
else
    [~,xmlBase,~] = fileparts(xmlParams.xmlFile);
    sqlDb.insert('experiment_info', {'start_time', 'task', 'session', 'animal', 'behavior_output_name'}, {datestr(now, 'yyyy-mm-dd HH:MM:SS'), xmlBase, sessionInfo, params.SubjectID, outfilename})
end

end