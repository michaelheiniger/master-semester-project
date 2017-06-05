function [timingEst, frequencyOffsetEst, initialPhaseEst] = timingAndFrequencyOffsetEstimation1(systemConfig, signal, stsTime, ltsTime, ca)
%TIMINGANDFREQUENCYOFFSETESTIMATION1 Find the beginning of the frame, the
%frequency offset and intial phase.
% Returned values are:
% - timingEst: estimate of the index of the first sample of the frame
% - frequencyOffsetEst: time-dependent frequency offset due different to
% carrier frequency at transmitter and receiver
% - intialPhaseEst: estimate intial phase offset (i.e. does not depend on time)

disp('Timing and Frequency Offset Estimation using STS/LTS and OFDM demodulation')

sc = systemConfig;

Ts = 1/sc.Fs;

% Coarse timing estimation
% Cross-correlation between 10 clean shorts STS with the signal
tenStsTime = repmat(stsTime,10,1); 
[R, ~] = xcorr(signal, repmat(stsTime,10,1));
R = R(length(signal):end);
R = R / sum(abs(tenStsTime).^2); % Normalization
[~, coarseTau] = max(abs(R))

plotSignalMagnitude(R, 'Offset [samples]', 'Cross-correlation (signal, 10 STS)');

% coarseTimingWindowLength = sc.CPLength+20; % maximum length of ISI 
coarseTimingWindowLength = sc.CPLength; % maximum length of ISI 

% Indices of samples candidate for first frame sample
fineTimingWindow = (coarseTau-coarseTimingWindowLength:coarseTau);
% fineTimingWindow = (coarseTau-coarseTimingWindowLength:coarseTau+coarseTimingWindowLength);

frameOffsetStartTwoLts = 10*length(stsTime) + sc.twoLtsCpLength;

allMSEs = zeros(1,length(fineTimingWindow));
allMAEs = zeros(1,length(fineTimingWindow));

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
    t1 = (0:length(ltsTimeRx1)-1)*Ts;
    t2 = (length(ltsTimeRx1):2*length(ltsTimeRx1)-1)*Ts;
    ltsTimeRxCorrected1 = ltsTimeRx1 .* exp(-1j*2*pi*coarseCFO(d)*t1).';
    ltsTimeRxCorrected2 = ltsTimeRx2 .* exp(-1j*2*pi*coarseCFO(d)*t2).';
    
    % Fine CFO estimation and correction
    fineCFO(d) = 1/(2*pi*length(ltsTime)*Ts)*angle(doInnerProduct(ltsTimeRxCorrected2, ltsTimeRxCorrected1));
    ltsTimeRxCorrected1 = ltsTimeRxCorrected1 .* exp(-1j*2*pi*fineCFO(d)*t1).';
    ltsTimeRxCorrected2 = ltsTimeRxCorrected2 .* exp(-1j*2*pi*fineCFO(d)*t2).';
    
    % Initial phase estimation and correction
    intialPhaseEst1 = angle(doInnerProduct(ltsTimeRxCorrected1, ltsTime));
    initialPhase(d) = intialPhaseEst1;
    intialPhaseEst2 = angle(doInnerProduct(ltsTimeRxCorrected2, ltsTime));
    ltsTimeRxCorrected1 = ltsTimeRxCorrected1 .* exp(-1j*intialPhaseEst1).';
    ltsTimeRxCorrected2 = ltsTimeRxCorrected2 .* exp(-1j*intialPhaseEst2).';
    
    % Go to frequency domain 
    ltsFreqRxCorrected1 = fftshift(fft(ltsTimeRxCorrected1,64));
    ltsFreqRxCorrected2 = fftshift(fft(ltsTimeRxCorrected2,64));
    ltsFreq = fftshift(fft(ltsTime,64));
    
    % Channel estimation over LTS1: it is protected from ISI due to
    % multipath by the cyclic prefix and from ISI due to late sync by LTS2
    lambdas = [ltsFreqRxCorrected1(7:32);ltsFreqRxCorrected1(34:59)] ./ [ltsFreq(7:32);ltsFreq(34:59)];

    zeroSubcarriersTop = ltsFreqRxCorrected1(1:6);
    zeroSubcarriersBottom = ltsFreqRxCorrected1(end-4);
    noiseVariance = 1/sc.numTotalCarriers*var([zeroSubcarriersTop; zeroSubcarriersBottom]);
    lambdasMMSE = channelEstimation([ltsFreqRxCorrected2(7:32);ltsFreqRxCorrected2(34:59)], [ltsFreqRxCorrected1(7:32);ltsFreqRxCorrected1(34:59)], [ltsFreq(7:32);ltsFreq(34:59)], sc.numUsedCarriers, noiseVariance);

    % Equalization of LTS2
    ltsFreqRxCorrected2(7:32) = ltsFreqRxCorrected2(7:32) ./ lambdasMMSE(1:26);
    ltsFreqRxCorrected2(34:59) = ltsFreqRxCorrected2(34:59) ./ lambdasMMSE(27:end);

    ltsFreqRxCorrected1(7:32) = ltsFreqRxCorrected1(7:32) ./ lambdasMMSE(1:26);
    ltsFreqRxCorrected1(34:59) = ltsFreqRxCorrected1(34:59) ./ lambdasMMSE(27:end);

    
    % Compute the error between the reference LTS and the demodulated LTS2,
    % Compare MSE vs MAE !
    allMSEs(d) = mean(abs([ltsFreq(7:32);ltsFreq(34:59); ltsFreq(7:32);ltsFreq(34:59)]-[ltsFreqRxCorrected1(7:32);ltsFreqRxCorrected1(34:59);ltsFreqRxCorrected2(7:32);ltsFreqRxCorrected2(34:59)]).^2); % Mean Squared Error 
%     allMSEs(d) = mean(abs(ltsFreq-ltsFreqRxCorrected2).^2); % Mean Squared Error 
%     allMAEs(d) = mean(abs(ltsFreq-ltsFreqRxCorrected2)); % Mean Absolute Error
    allMAEs(d) = mean(abs([ltsFreq(7:32);ltsFreq(34:59)]-[ltsFreqRxCorrected2(7:32);ltsFreqRxCorrected2(34:59)])); % Mean Absolute Error
end

threshold = 0.25;
quantileError = quantile(allMSEs, [threshold, 1])
% Take the sample which has an error below the treshold and which has the highest index 
posQuantile = find(allMSEs <= quantileError(1), 1, 'last')
fineTauQuantile = fineTimingWindow(posQuantile)

figure;
hax=axes;
stem(fineTimingWindow, abs(allMSEs),'.')
HL1 = quantileError(1);
line(get(hax,'XLim'),[HL1 HL1], 'Color', [1 0 0]);
title('Total error on LTS vs offset');

[minMSE, posMinMSE] = min(allMSEs);
[minMAE, posMinMAE] = min(allMAEs);
disp(['Minimum MSE: ', num2str(minMSE)]);
disp(['Minimum MAE: ', num2str(minMAE)]);

fineTauMinMSE = fineTimingWindow(posMinMSE);
fineTauMinMAE = fineTimingWindow(posMinMAE);
fineTau = fineTauQuantile;
% fineTau = fineTimingWindow(posMinMSE)

position = posQuantile

% Begining of the frame: take into account the CA code
frameBeginning = fineTau-length(ca);
disp(['Frame beginning est:', num2str(frameBeginning)]);

% Timing estimates is the estimate of the first sample of the frame
timingEst = frameBeginning;

% Total frequency offset is the sum of the coarse and fine estimations
frequencyOffsetEst = coarseCFO(position) + fineCFO(position);
disp(['Coarse CFO: ', num2str(coarseCFO(position))]);
disp(['Fine CFO: ', num2str(fineCFO(position))]);
disp(['Total CFO: ', num2str(frequencyOffsetEst)]);

% Initial phase offset wrt the beginning of the frame
% Note: Initial phase offset is computed wrt LTS1 which is preceded by
% CA, 10*STSs and the CP of LTS1 and LTS2
offsetToFrameBeginning = length(ca)+length(tenStsTime)+sc.twoLtsCpLength;
initialPhaseEst = angle(exp(1j*initialPhase(position)-1j*(2*pi*frequencyOffsetEst*offsetToFrameBeginning*Ts)));
disp(['Initial phase:', num2str(initialPhaseEst)]);

end

