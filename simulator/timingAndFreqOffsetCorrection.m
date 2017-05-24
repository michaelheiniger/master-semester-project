function [frameCorrected, frequencyOffsetEst] = timingAndFreqOffsetCorrection(signal, stsTime, ltsTime, numTotalCarriers, CPLength, twoLtsCpLength, numOFDMSymbolsPerFrame, Fs)
%TIMINGANDFREQOFFSET Summary of this function goes here
%   Detailed explanation goes here

global decoder;

frequencyOffsetEst = 0;

Ts = 1/Fs;

lengthStsTime = length(stsTime);
lengthLtsTime = length(ltsTime);

% Coarse timing estimation
% Cross-correlation between 10 clean shorts STS with the signal
tenStsTime = repmat(stsTime,10,1); 
[R, ~] = xcorr(signal, repmat(stsTime,10,1));
R = R(length(signal):end);
R = R / sum(abs(tenStsTime).^2); % Normalization
[~, coarseTau] = max(abs(R))

figure;
plot(abs(R));
xlabel('Offset [samples]');
ylabel('Magnitude');
title('Cross-correlation (signal, 10 STS)');


coarseTimingWindowLength = CPLength; % maximum length of ISI 

% Indices of samples candidate for first frame sample
fineTimingWindow = (coarseTau-coarseTimingWindowLength:coarseTau);
frameOffsetStartTwoLts = 10*length(stsTime) + twoLtsCpLength;

totalErrors = zeros(1,length(fineTimingWindow));

% Store the CFO estimations for every tentative timing
coarseCFO = zeros(1,length(fineTimingWindow));
fineCFO = zeros(1,length(fineTimingWindow));
initialPhase = zeros(1,length(fineTimingWindow));

% Loop over all timing candidates
for d = 1:length(fineTimingWindow)
    ltsTimeRx1 = signal(fineTimingWindow(d)+frameOffsetStartTwoLts:fineTimingWindow(d)+frameOffsetStartTwoLts+length(ltsTime)-1); % LTS 1
    ltsTimeRx2 = signal(fineTimingWindow(d)+frameOffsetStartTwoLts+length(ltsTime):fineTimingWindow(d)+frameOffsetStartTwoLts+2*length(ltsTime)-1); % LTS 2
        
    stsTime8 = signal(fineTimingWindow(d)+7*length(stsTime):fineTimingWindow(d)+8*length(stsTime)-1);
    stsTime9 = signal(fineTimingWindow(d)+8*length(stsTime):fineTimingWindow(d)+9*length(stsTime)-1);
    
    % Coarse CFO estimation and correction
    coarseCFO(d) = 1/(2*pi*length(stsTime)*Ts)*angle(doInnerProduct(stsTime9, stsTime8));
    t = (0:length(ltsTimeRx1)-1)*Ts;
    ltsTimeRxCorrected1 = ltsTimeRx1 .* exp(-1j*2*pi*coarseCFO(d)*t).';
    ltsTimeRxCorrected2 = ltsTimeRx2 .* exp(-1j*2*pi*coarseCFO(d)*t).';
    
    % Fine CFO estimation and correct
    fineCFO(d) = 1/(2*pi*length(ltsTime)*Ts)*angle(doInnerProduct(ltsTimeRxCorrected2, ltsTimeRxCorrected1));
    t = (0:length(ltsTimeRxCorrected1)-1)*Ts;
    ltsTimeRxCorrected1 = ltsTimeRxCorrected1 .* exp(-1j*2*pi*fineCFO(d)*t).';
    ltsTimeRxCorrected2 = ltsTimeRxCorrected2 .* exp(-1j*2*pi*fineCFO(d)*t).';
    
    % Find and correct initial estimate
    intialPhaseEst1 = angle(doInnerProduct(ltsTimeRxCorrected1, ltsTime));
    initialPhase(d) = intialPhaseEst1;
    intialPhaseEst2 = angle(doInnerProduct(ltsTimeRxCorrected2, ltsTime));
    ltsTimeRxCorrected1 = ltsTimeRxCorrected1 .* exp(-1j*intialPhaseEst1).';
    ltsTimeRxCorrected2 = ltsTimeRxCorrected2 .* exp(-1j*intialPhaseEst2).';
    
    ltsFreqRxCorrected1 = fftshift(fft(ltsTimeRxCorrected1,64));
    ltsFreqRxCorrected2 = fftshift(fft(ltsTimeRxCorrected2,64));
    ltsFreq = fftshift(fft(ltsTime,64));
    
    % Channel estimation over LTS1: it is protected from ISI due to
    % multipath by the cyclic prefix and from ISI due to late sync by LTS2
    lambdas = [ltsFreqRxCorrected1(7:32);ltsFreqRxCorrected1(34:59)] ./ [ltsFreq(7:32);ltsFreq(34:59)];

    % Equalization of LTS2
    ltsFreqRxCorrected2(7:32) = ltsFreqRxCorrected2(7:32) ./ lambdas(1:26);
    ltsFreqRxCorrected2(34:59) = ltsFreqRxCorrected2(34:59) ./ lambdas(27:end);

%     innerProducts(d) = doInnerProduct(ltsFreqRxCorrected2,ltsFreqRxCorrected2);
    
    % Compute the error between the reference LTS and the demodulated LTS2
    totalErrors(d) = sum(abs(ltsFreq-ltsFreqRxCorrected2));
end

threshold = 0.25;
quantileError = quantile(totalErrors, [threshold, 1])
% Take the sample which has an error below the treshold and which has the highest index 
[~, posQuantile] = find(totalErrors <= quantileError(1), 1, 'last');
fineTauQuantile = fineTimingWindow(posQuantile);

figure;
stem(fineTimingWindow, abs(totalErrors))
title('Total error on LTS vs offset');

    
[~, posMinError] = min(totalErrors);
fineTauMinError = fineTimingWindow(posMinError)
fineTau = fineTauQuantile
position = posQuantile;

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
    [intialPhaseEst, frequencyOffsetEst] = frequencyOffsetEstimation(frame, stsTime, ltsTime, Fs);
    
    totalCFO = coarseCFO(position) + fineCFO(position);
    
    % Correct frequency offset
    t = (0:length(frame(1+10*lengthStsTime+twoLtsCpLength:end))-1)*Ts;
    frameCorrected = [frame(1:10*lengthStsTime+twoLtsCpLength); frame(1+10*lengthStsTime+twoLtsCpLength:end) .* exp(-1j*(2*pi*totalCFO*t+initialPhase(position))).'];
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Old code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lengthStsTime = length(stsTime);
% lengthLtsTime = length(ltsTime);
% 
% if decoder.upsample
%     signal = resample(signal, decoder.USF, 1);
%     % Timing estimation - find the begining of the next received OFDM frame
%     fineTau = timingSync(signal, resample(stsTime, decoder.USF, 1), resample(ltsTime, decoder.USF, 1), twoLtsCpLength, CPLength);
% else
%     % Timing estimation - find the begining of the next received OFDM frame
%     fineTau = timingSync(signal, stsTime, ltsTime, twoLtsCpLength, CPLength);    
% end
% 
% % Extract current frame
% beginingFrame = fineTau;
% signalCut = signal(beginingFrame:end);
% if decoder.upsample
%     signalCut = downsample(signalCut, decoder.USF);
% end
% lengthFrame = 10*lengthStsTime + twoLtsCpLength + 2*lengthLtsTime + numOFDMSymbolsPerFrame*(CPLength+numTotalCarriers);
% frame = signalCut(1:lengthFrame);
% 
% figure;
% hax=axes;
% plot(abs(signal))
% title('Abs val of received signal');
% VL1 = fineTau;
% line([VL1 VL1],get(hax,'YLim'), 'Color', [0 1 0]);
% VL2 = fineTau+lengthFrame-1;
% line([VL2 VL2],get(hax,'YLim'), 'Color', [0 1 0]);
% 
% figure;
% plot(abs(frame));
% title('Received frame');
% 
% if decoder.cfoCorrection
%     [intialPhaseEst, frequencyOffsetEst] = frequencyOffsetEstimation(frame, stsTime, ltsTime, Fs);
% 
% %     frequencyOffsetEst = 237 %TMP: Used to set manually the residual CFO error
%     
%     % Correct frequency offset
%     t = (0:length(frame(1+10*lengthStsTime+twoLtsCpLength:end))-1)*Ts;
%     frameCorrected = [frame(1:10*lengthStsTime+twoLtsCpLength); frame(1+10*lengthStsTime+twoLtsCpLength:end) .* exp(-1j*(2*pi*frequencyOffsetEst*t+intialPhaseEst)).'];
% else
%     frameCorrected = frame;
% end
% 
% if decoder.upsample
%     figure;
%     hax=axes;
%     plot(abs(signal))
%     title('Abs val of upsampled received signal');
%     VL1 = fineTau;
%     line([VL1 VL1],get(hax,'YLim'), 'Color', [0 1 0]);
%     VL2 = fineTau+decoder.USF*lengthFrame-1;
%     line([VL2 VL2],get(hax,'YLim'), 'Color', [0 1 0]);
% end

% figure;
% plot(abs(frameCorrected));
% title('Corrected received frame');

end