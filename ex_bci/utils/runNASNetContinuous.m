function waveformPass = runNASNetContinuous(w1, b1, w2, b2, wvForms, gamma)

gammaTransformed = log(gamma/(1-gamma));
nasTransformedValue = max(0,wvForms*w1+b1')*w2 + b2;

waveformPass = nasTransformedValue>gammaTransformed;