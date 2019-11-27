function n = pullOutNumber(s)

s = s(3:strfind(s,'cd/m2')-1);

switch s(end)
    case ' '
        n = str2double(s);
    case 'm'
        n = str2double(s(1:end-1))/1000;
    case 'k'
        n = str2double(s(1:end-1))*1000;
    otherwise
        disp('Problem');
end

