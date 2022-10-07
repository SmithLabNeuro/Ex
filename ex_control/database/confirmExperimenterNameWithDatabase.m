function confirmedExperimenter = confirmExperimenterNameWithDatabase(experimenter)
% adds a new subject to the database--only needs to run once per subject!
% (or... should...)
global sqlDb
if isempty(sqlDb)
    error('No database linked, so nothing could be read from/written to database.')
end 

experimenterCheck = sqlDb.fetch(sprintf('SELECT experimenter FROM experiment_session WHERE experimenter=''%s''', experimenter));
if isempty(experimenterCheck)
    addCheck = questdlg(sprintf('%s is a new main experimenter. Is this correct?', experimenter));
    switch addCheck
        case 'Yes'
            confirmedExperimenter = experimenter;
        case {'No', 'Cancel'}
            correctName = inputdlg('Enter correct experimenter:');
            if isempty(correctName)
                error('Experimenter needed!')
            else
                correctName = correctName{1};
            end
            correctName = regexprep(lower(correctName),'\s',''); %enforce lowercase / no whitespace
            confirmedExperimenter = confirmExperimenterNameWithDatabase(correctName);
        otherwise
            error('Experimenter needed!')
    end
else    
    confirmedExperimenter = experimenter;
end
end