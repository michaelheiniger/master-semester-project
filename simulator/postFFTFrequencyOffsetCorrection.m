function [ofdmSymbolsCorrected] = postFFTFrequencyOffsetCorrection(ofdmSymbols, channelCoefficients, pilots1, pilots2, numDataCarriers, numZeros)
%POSTFFTFREQUENCYOFFSETESTIMATION Summary of this function goes here
%   Detailed explanation goes here

pilots1Rx = ofdmSymbols(1,2:end);
pilots2Rx = ofdmSymbols(end,2:end);

phase1 = angle(pilots1Rx.*conj(pilots1));
phase2 = angle(pilots2Rx.*conj(pilots2));

slopes = abs(phase1-phase2)/numDataCarriers;

figure;
plot(slopes,'.');

K = repmat((-numDataCarriers/2:numDataCarriers/2-1).', 1, size(pilots1Rx,2));

% size(K)
% size(repmat(slopes, numDataCarriers, 1))
% size(ofdmSymbols(:,2:end))

ofdmSymbolsCorrected = ofdmSymbols(:,2:end) .* exp(1j*K.*slopes);


end

