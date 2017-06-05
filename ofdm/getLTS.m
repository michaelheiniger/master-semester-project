function [ltsTime, ltsFreq] = getLTS()
%GETLTS IEEE 802.11a Long Training Sequence
%  returns the LTS in time and frequency domain

NFFT = 64;

% Long Training Sequence (LTS) in frequency domain
% it is already extended to 64 and scaled as in 802.11a
load('training-sequences/ltsFreq.mat');

% LTS in time domain
ltsTime = ifft(fftshift(ltsFreq), NFFT);

end

