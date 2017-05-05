function [frameCorrected, frequencyOffsetEst] = timingAndFreqOffsetCorrection(signal, stsTime, ltsTime, numTotalCarriers, CPLength, twoLtsCpLength, numOFDMSymbolsPerFrame)
%TIMINGANDFREQOFFSET Summary of this function goes here
%   Detailed explanation goes here

global decoder;

frequencyOffsetEst = 0;

Fs = 0.5e6;
Ts = 1/Fs;

lengthStsTime = length(stsTime);
lengthLtsTime = length(ltsTime);

if decoder.upsample
    signal = resample(signal, decoder.USF, 1);
    % Timing estimation - find the begining of the next received OFDM frame
    fineTau = timingSync(signal, resample(stsTime, decoder.USF, 1), resample(ltsTime, decoder.USF, 1), twoLtsCpLength);
else
    % Timing estimation - find the begining of the next received OFDM frame
    fineTau = timingSync(signal, stsTime, ltsTime, twoLtsCpLength);    
end

% Extract current frame
beginingFrame = fineTau;
signalCut = signal(beginingFrame:end);
if decoder.upsample
    signalCut = downsample(signalCut, decoder.USF);
end
lengthFrame = 10*lengthStsTime + twoLtsCpLength + 2*lengthLtsTime + numOFDMSymbolsPerFrame*(CPLength+numTotalCarriers);
frame = signalCut(1:lengthFrame);

figure;
hax=axes;
plot(abs(signal))
title('Abs val of received signal');
VL1 = fineTau;
line([VL1 VL1],get(hax,'YLim'), 'Color', [0 1 0]);
VL2 = fineTau+lengthFrame-1;
line([VL2 VL2],get(hax,'YLim'), 'Color', [0 1 0]);

figure;
plot(abs(frame));
title('Received frame');

if decoder.cfoCorrection
    [intialPhaseEst, frequencyOffsetEst] = frequencyOffsetEstimation(frame, stsTime, ltsTime);

%     frequencyOffsetEst = 237 %TMP: Used to set manually the residual CFO error
    
    % Correct frequency offset
    t = (0:length(frame(1+10*lengthStsTime+twoLtsCpLength:end))-1)*Ts;
    frameCorrected = [frame(1:10*lengthStsTime+twoLtsCpLength); frame(1+10*lengthStsTime+twoLtsCpLength:end) .* exp(-1j*(2*pi*frequencyOffsetEst*t+intialPhaseEst)).'];
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
plot(abs(frameCorrected));
title('Corrected received frame');

end