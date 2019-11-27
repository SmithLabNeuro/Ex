function mov=loadcx_movie(filename)

fid=fopen(filename, 'r','l');
notes=setstr(fread(fid, 10, 'char'))'
dmns=fread(fid, 4, 'uint16')';

if dmns(4) ~= 1
  dmns(4)=dmns(4)+1; 
end

%dmns
dmns

mov = zeros(dmns(3),dmns(2),dmns(4),'uint8');

for i=1:dmns(4)
    imgmtx = fread(fid, dmns(2)*dmns(3), 'uchar');
	imgmtx = reshape(imgmtx, dmns(2), dmns(3));
%    i
notes = setstr(fread(fid,10,'char'));
dmns = fread(fid,4,'uint16');
    mov(:,:,i) = uint8(imgmtx)';
%    mov(:,:,i) = uint8(imgmtx);
end

fclose(fid);
