
%waveforms = randn(1000,52);
clear all
waveforms = [{}];
for n = 1:200
    waveforms{n} = randn(5,52);
end
W1 = randn(52,52);
B1 = randn(52,1024); %equal to buffer size
W2 = randn(52,2);
B2 = randn(2,1);
bound = B2(1)-B2(2);
zerovec = zeros(52,5);
relu = @(x)(max(zerovec,x));
a=tic;
goodcounts = zeros(200,1);
for n = 1:length(goodcounts)
    thiswaves = waveforms{n};
    numwaves = size(thiswaves,1);
    temp = W2'*max(0,W1*thiswaves'+B1(:,1:numwaves));
    goodcounts(n)=sum((temp(2,:)-temp(1,:))>bound);%this approximates sum((W2'*max(0,W1*waveforms{n}'+b1))>4);%
end%(n)=sum((W2'*max(zerovec,W1*waveforms{n}'+b1)));%

toc(a)