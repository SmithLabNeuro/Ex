frames = 25;
dim = 400;

for i = 1:60
  mov = gaussNoiseMtx(frames,8,dim,1);
  
  save(sprintf('noise%i',i),'mov');
end