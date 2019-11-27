function savecx_movie(filename, notes, dmns, images)

% SAVECX_MOVIE(filename, notes, dmns, imgmtx)
%       save the image files as a cortex readable image file.
%       filename,       path should be included
%       notes,          maximume 10 characters
%       dmns=[depth, x, y, nframes], in which
%               depth,          bitmap depth (1,2,4, or 8)
%               x,              x dimension of the image
%               y,              y dimension of the image
%               nframes,        number of frames in the movie (1's based)
%       imgmtx,    A cell array of images, will be rounded to the range of 0-255.
%
% By Yi-Xiong Zhou on 4-9-96
% Modified by Brian Potetz to save CTX movies
%

% Truncate/extend the "notes" field to 10 characters
if size(notes,2)>10
    notes = notes([1:10]);
else
    notes = [notes, zeros(1,10-size(notes,2))];
end

% Make sure the dimensions line up
x=dmns(2); y=dmns(3); nf=dmns(4);
imx = size(images{1},2);
imy = size(images{1},1);
if (nf ~= length(images))
    disp(sprintf('ERROR: Number of frames (%d) does not match length of cell array (%d)', nf, length(images)));
    return;
end

% Check all image sizes
for i=2:nf
    if (size(images{i},2) ~= imx) | (size(images{i},1) ~= imy)
        disp(sprintf('ERROR: Size of image #%d does not match size of image #1', i));
        return;
    end
end

% Cortex counts frames starting at zero.
% This routine expects the user to count the frames starting at one, so
%  we translate the number of frames to zero-based here.
dmns(4) = max(0, dmns(4)-1);


fid=fopen(filename, 'w');

[fn, pp, ar]=fopen(fid);
if strcmp(ar, 'ieee-be')
	tmp=floor(dmns/256);
	dmns=(dmns-tmp*256)*256+tmp;
elseif ~(strcmp(ar, 'ieee-le') | strcmp(ar,'ieee-le.l64'))
	disp('unknow file format. find out the byte switch requirement!')
	disp('use unswitch as default') 
end

for i=1:nf
	fwrite(fid, notes, 'char');
	fwrite(fid, dmns, 'uint16');
	if i==1, dmns(4)=0; end
% 	fwrite(fid, imgmtx(:, 1+(i-1)*y:i*y), 'uchar');
	fwrite(fid, images{i}', 'uchar');
%     nshow(images{i})
end
fclose(fid);




