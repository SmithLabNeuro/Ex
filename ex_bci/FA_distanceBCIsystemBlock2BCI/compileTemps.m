clear all
path = '../../../../../Documents/bcitemp/';
savepath = '../../../../../bciData/pepelepew/';
filename = 'Pe190307_test';
file=dir(path);

goodfiles = [{''}];
filenums = [];
count = 1;
for n = 1:length(file)
    if ~isempty(strfind(file(n).name,'onlinedat'))
        goodfiles{count} = file(n).name;
        [~,name,ext]=fileparts(file(n).name((length('onlinedat')+1):end));
        filenums(count) = str2double(name);
        count = count+1;
    end
end

[~,inds] = sort(filenums,'ascend');
orderedfiles = goodfiles(inds);

for n = 1:length(orderedfiles)
    load([path,orderedfiles{n}])
    onlinedat(n) = temponlinedat;
end


testfile = [savepath,filename,'onlinedat1.mat'];
postind = 1;
while exist(testfile,'file')
postind = postind+1;
testfile = [savepath,filename,'onlinedat',num2str(postind),'.mat'];
end
save(testfile,'onlinedat')