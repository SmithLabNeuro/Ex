load Pe190115_s514a_distanceStabilityBCIcalib_0001_lowdHmm0onlinedat1
nosham= onlinedat([onlinedat.shamflag]==0);
nosham = nosham(61:end);
lens = [];for n =1:length(nosham);lens(n)=length(nosham(n).trackdist);end
shamdat2 = nosham(lens==76);
shamdat2 = shamdat2([shamdat2.correctflag] == 1);
shamdat = [shamdat shamdat2]