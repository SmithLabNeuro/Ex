function frame = ctxNoise(nframes,pixelsize,len,centerLum,lumChange,filename)
%function frame = ctxNoise(nframes,pixelsize,len,centerLum,lumChange,filename)
%    nframes: number of frames in the ctx movie
%    pixelsize: width of the noise blocks, in pixels
%    len: width of the entire movie, in pixels
%    centerLum: the mean luminance value 
%    lumChange: the luminance values will range by this amount 
%               (centerLum +- lumChange)
%    filename: the filename of the output ctx file (e.g. test.ctx)
%
%    Example: generate a movie that's 200x200 with 40 pixel blocks, 1000
%    frames which has values from 98-138, and output it to file 'blah.ctx'
%       >> c = ctxNoise(1000,40,200,118,20,'blah.ctx');
%
%    To view the frames outputted by c = ctxNoise(...)
%       >> colormap gray % do this once
%       >> imagesc(c{57},[128 255]) % or whatever range you want ([0 255])
%
%    To load the ctx file later in matlab:
%       >> mov = loadcx_movie('test.ctx');
%       >> imagesc(mov(:,:,57),[128 255]);

  rep = 1;

    frame = zeros(len,len,nframes);

    len = len / pixelsize;
        
    for i = 1:nframes
        pixels = rand(len);
        thisframe = ...
            reshape(...
                    permute(...
                            repmat(...
                                   reshape(...
                                           repmat(...
                                                  reshape(pixels,1,len*len),...
                                                  [1,pixelsize,1] ...
                                                  ), ...
                                           size(pixels,1), ...
                                           len*pixelsize ...
                                           )', ...
                                    [1 1 pixelsize] ...
                                    ),...
                             [3 1 2] ...
                             ), ...
                     len*pixelsize, ...
                     len*pixelsize ...
                     );
        frame(:,:,i) = circshift(thisframe,floor(rand(2,1)*pixelsize));
    end
    
    m = frame * lumChange * 2 - lumChange;
    m = m + centerLum;

    clear frame
    
    frame = cell(nframes*rep,1);
    
    for i = 1:nframes
        for j = 1:rep
            f = uint8(m(:,:,i));
            frame{(i-1)*rep+j} = f;            
        end
    end

    savecx_movie(filename,filename,[8 len*pixelsize len*pixelsize nframes],frame);