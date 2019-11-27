function out = driftchoiceextractparam(dat,str,pickwhich)
pickwhichflag = exist('pickwhich','var');
out = nan(length(dat),1);
for n = 1:length(dat)
    codes = dat(n).events;
    textcodes = codes(codes(:,2)>255,2)-256;
    text = char(textcodes)';
    startstrind= strfind(text,str)+length(str);
    if pickwhichflag == 1
        startstrind = startstrind(min(pickwhich,length(startstrind)));
    end
    startstrindcol = strfind(text(startstrind(1):end),';')-1;
    out(n) = str2double(strsplit(text(startstrind(1):startstrind(1)+startstrindcol-1)));
end
end