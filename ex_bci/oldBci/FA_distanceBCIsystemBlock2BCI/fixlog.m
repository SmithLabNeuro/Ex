%%fix this nev

filename='../../../../../bciData/pepelepew/Pe181130_s488a_distanceStabilityBCIcalib_0001_lowdHmm0_log1.txt';
fileid = fopen(filename);
temp = 0;
a=[{}];
count = 1;
while temp ~=-1
temp=fgetl(fileid);
a{count} = temp;
count = count+1;
end
fclose(fileid)
a2=[{}];
shamnum = '0';
fid = fopen('myfix.txt','w')
for n = 1:length(a)
    n
    asplit = strsplit(a{n},'\t');
    if n==1
        %a2{n} = asplit;
        fprintf(fid,[a{1},'\n']);
    elseif length(asplit)==7
        if n ~= length(a)
            anext = strsplit(a{n+1},'\t');
            numstart = strfind(asplit{2},anext{1});
            if length(numstart)>1
                numstart = numstart(end);
            end
            trialnum = asplit{2}(numstart(1):end);
            shamnum = asplit{2}(1:(numstart(1)-1));
            if strcmp(shamnum,'-1')
                shamnum = '0';
            end
            asplit(1) = [];
            asplit{1} = trialnum;
            asplit{end+1}=shamnum;
            %str = sprintf([asplit{1},'\t',asplit{2},'\t',asplit{3},'\t',asplit{4},'\t',asplit{5},'\t',asplit{6},'\t',asplit{7}]);
            %a2{n} = str;
            fprintf(fid,[asplit{1},'\t',asplit{2},'\t',asplit{3},'\t',asplit{4},'\t',asplit{5},'\t',asplit{6},'\t',asplit{7},'\n']);
        end
    else
        fprintf(fid,[a{n} ,'\t',shamnum,'\n']);
        %str = sprintf([a{n} ,'\t',shamnum]);
        %a2{n} = str;
    end
end
fclose(fid)

