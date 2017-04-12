function [frameCorrected, frequencyOffsetEst] = timingAndFreqOffsetCorrection(signal, shortTrainingSeqIfft, longTrainingSeqIfft, numTotalCarriers, CPLength, preamble2CPLength, numOFDMSymbolsPerFrame)
%TIMINGANDFREQOFFSET Summary of this function goes here
%   Detailed explanation goes here

frequencyOffsetEst = 0;

% receivedSignal = signal;

global decoder;

Fs = 0.5e6;
Ts = 1/Fs;

lengthSTSIfft = length(shortTrainingSeqIfft);
lengthLTSIfft = length(longTrainingSeqIfft);

if decoder.upsample
    signal = resample(signal, decoder.USF, 1);
    % Timing estimation - find the begining of the next received OFDM frame
    fineTau = timingEstimation(signal, resample(shortTrainingSeqIfft, decoder.USF, 1));
else
    % Timing estimation - find the begining of the next received OFDM frame
    fineTau = timingEstimation(signal, shortTrainingSeqIfft);    
end

% Extract current frame
beginingFrame = fineTau;
signalCut = signal(beginingFrame:end);
if decoder.upsample
    signalCut = downsample(signalCut, decoder.USF);
end
lengthFrame = 10*lengthSTSIfft + preamble2CPLength + 2*lengthLTSIfft + numOFDMSymbolsPerFrame*(CPLength+numTotalCarriers);
frame = signalCut(1:lengthFrame);

if decoder.frequencyCorrection
    [intialPhaseEst, frequencyOffsetEst] = frequencyOffsetEstimation(frame, shortTrainingSeqIfft, longTrainingSeqIfft);

    % Correct frequency offset
    t = (0:length(frame(1+10*lengthSTSIfft+preamble2CPLength:end))-1)*Ts;
    frameCorrected = [frame(1:10*lengthSTSIfft+preamble2CPLength); frame(1+10*lengthSTSIfft+preamble2CPLength:end) .* exp(-1j*(2*pi*frequencyOffsetEst*t+intialPhaseEst)).'];
else
    frameCorrected = frame;
end

if decoder.upsample
    figure;
    hax=axes;
    plot(abs(signal))
    title('Abs val of upsampled received signal');
    VL1 = fineTau;
    line([VL1 VL1],get(hax,'YLim'), 'Color', [0 1 0]);
    VL2 = fineTau+decoder.USF*lengthFrame-1;
    line([VL2 VL2],get(hax,'YLim'), 'Color', [0 1 0]);
end

figure;
% hax=axes;
plot(abs(frameCorrected)); %/sum(abs(normalization(1:length(R))).^2));
title('Corrected received signal');

end