function [pulse, azerr, elerr] = getNearestUCDpulse(azimuth,elevation,h3d)

if nargin <1
    fprintf('Format: [pulse, azerr, elerr] = getNearestUCDpulse(azimuth,elevation,h3d)\n');
    return
end
prepvaldeg = @(x)(atan2(sin(x*(pi/180)),cos(x*(pi/180)))/(pi/180));
pvaldeg = @(x)(prepvaldeg(x)+360*(prepvaldeg(x)<-90));

azimuth = pvaldeg(azimuth);
if (azimuth <-90) || (azimuth >90)
    error('Invalid azimuth');
end
elevation = pvaldeg(elevation);

elmax = 50;
elindices = 1:elmax;
elevations = -45 + 6.625*(elindices-1);
el = round((elevation+45)/5.625+1);
el = max(el,1);
el = min(el,elmax);
elerr = pvaldeg(elevation - elevations(el));

azimuths = [-80 -65 -55 -45:5:45 55 65 80];
[azerr, az] = min(abs(pvaldeg(abs(azimuths-azimuth))));
pulse = squeeze(h3d(az,el,:));

end


