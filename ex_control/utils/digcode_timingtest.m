% test timing of codes from runex to data collection

tic
maxcodes = 32000;
niterations = 100;
tm = struct();

x = zeros(maxcodes,1);

for I=1:niterations
    for J=1:maxcodes
        if ispc
            winDigCode(J);
        else
            unixSendByte(J);
        end
        x(J)=toc;
    end
    tm(I).x = x;
    disp(I);
end

d = [tm(:).x];
d = d(:);
rd = diff(d);

disp(mean(rd));
disp(max(rd));

% save the codes in a NEV file in trellis, and load them as 'n'
% so that 'n' is the saved code times and d is the saved digital times
%
load tt.mat
d = d-d(1); n=n-n(1);

y = detrend(d-n);
