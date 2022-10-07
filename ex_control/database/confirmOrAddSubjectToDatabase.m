function confirmedName = confirmOrAddSubjectToDatabase(subjectName)
% adds a new subject to the database--only needs to run once per subject!
% (or... should...)
global sqlDb
if isempty(sqlDb)
    error('No database linked, so nothing could be read from/written to database.')
end 

confirmedName = subjectName;
subjectCheck = sqlDb.fetch(sprintf('SELECT name FROM animal WHERE name=''%s''', subjectName));
if isempty(subjectCheck)
    addCheck = questdlg(sprintf('%s not in database. Add?', subjectName));
    switch addCheck
        case 'Yes'
            birthYearStr = inputdlg(sprintf('What is %s''s birth year?', subjectName));
            birthYear = str2double(birthYearStr);
            sqlDb.insert('animal', {'name', 'birth_year'}, {subjectName, birthYear});
        case {'No', 'Cancel'}
            correctName = inputdlg('Enter correct name if desired:');
            if isempty(correctName)
                error('Database unusable with nonexistent animal')
            else
                correctName = correctName{1};
            end
            correctName = regexprep(lower(correctName),'\s',''); %enforce lowercase / no whitespace
            confirmedName = confirmOrAddSubjectToDatabase(correctName);
        otherwise
            error('Database unusable with nonexistent animal')
    end
end
end