function [ stsTime, stsFreq ] = getSTS()
%GETSTS IEEE 802.11a Short Training Sequence
%  returns the STS in time and frequency domain

NFFT = 64;

% Short Training Sequence (STS) in frequency domain
% it is already extended to 64 and scaled as in 802.11a
load('training-sequences/stsFreq.mat');

% One FFT period (N = 64 samples) contains 4 STSs
fourStsTime = ifft(fftshift(stsFreq), NFFT);

% STS in time domainm
stsTime = fourStsTime(1:16);

end

