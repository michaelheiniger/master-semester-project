function [tau] = timingSync(signal, stsTime, ltsTime, twoLtsCpLength)
%TIMINGESTIMATION provide timing estimate of the begining of the next
%received frame

global decoder;

Ts = 1/(0.5e6);

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


coarseWindowLength = 16; % length of ISI 
frameOffsetStartTwoLts = 10*length(stsTime) + twoLtsCpLength;
% fineWindow = (coarseTau-coarseWindowLength:coarseTau+coarseWindowLength);
fineWindow = (coarseTau-coarseWindowLength:coarseTau);
totalErrors = zeros(1,length(fineWindow));
innerProducts = zeros(1,length(fineWindow));
% normalization = sum(abs(twoLtsTime).^2);
for d = 1:length(fineWindow)
    ltsTimeRx1 = signal(fineWindow(d)+frameOffsetStartTwoLts:fineWindow(d)+frameOffsetStartTwoLts+length(ltsTime)-1); % LTS 1
    ltsTimeRx2 = signal(fineWindow(d)+frameOffsetStartTwoLts+length(ltsTime):fineWindow(d)+frameOffsetStartTwoLts+2*length(ltsTime)-1); % LTS 2
        
    stsTime8 = signal(fineWindow(d)+7*length(stsTime):fineWindow(d)+8*length(stsTime)-1);
    stsTime9 = signal(fineWindow(d)+8*length(stsTime):fineWindow(d)+9*length(stsTime)-1);
    
    % Coarse CFO estimation and correction
    coarseCFO = 1/(2*pi*length(stsTime)*Ts)*angle(doInnerProduct(stsTime9, stsTime8));
    t = (0:length(ltsTimeRx1)-1)*Ts;
    ltsTimeRxCorrected1 = ltsTimeRx1 .* exp(-1j*2*pi*coarseCFO*t).';
    ltsTimeRxCorrected2 = ltsTimeRx2 .* exp(-1j*2*pi*coarseCFO*t).';
    
    % Fine CFO estimation and correct
    fineCFO = 1/(2*pi*length(ltsTime)*Ts)*angle(doInnerProduct(ltsTimeRxCorrected2, ltsTimeRxCorrected1));
    t = (0:length(ltsTimeRxCorrected1)-1)*Ts;
    ltsTimeRxCorrected1 = ltsTimeRxCorrected1 .* exp(-1j*2*pi*fineCFO*t).';
    ltsTimeRxCorrected2 = ltsTimeRxCorrected2 .* exp(-1j*2*pi*fineCFO*t).';
    
    % Find and correct initial estimate
    intialPhaseEst1 = angle(doInnerProduct(ltsTimeRxCorrected1, ltsTime));
    intialPhaseEst2 = angle(doInnerProduct(ltsTimeRxCorrected2, ltsTime));
    ltsTimeRxCorrected1 = ltsTimeRxCorrected1 .* exp(-1j*intialPhaseEst1).';
    ltsTimeRxCorrected2 = ltsTimeRxCorrected2 .* exp(-1j*intialPhaseEst2).';
    
    ltsFreqRxCorrected1 = fftshift(fft(ltsTimeRxCorrected1,64));
    ltsFreqRxCorrected2 = fftshift(fft(ltsTimeRxCorrected2,64));
    ltsFreq = fftshift(fft(ltsTime,64));
    
    
    % Channel estimation over the first LTS: LTS1 is protected of ISI (due to multipath) by CP
    % and of ISI (due to late sync) by LTS2
    lambdas = [ltsFreqRxCorrected1(7:32);ltsFreqRxCorrected1(34:59)] ./ [ltsFreq(7:32);ltsFreq(34:59)];
%     lambdas = [ltsFreqRxCorrected2(7:32);ltsFreqRxCorrected2(34:59)] ./ [ltsFreq(7:32);ltsFreq(34:59)];
%     lambdas = channelEstimation([ltsFreqRxCorrected1, ltsFreqRxCorrected2]);
    ltsFreqRxCorrected2(7:32) = ltsFreqRxCorrected2(7:32) ./ lambdas(1:26);
    ltsFreqRxCorrected2(34:59) = ltsFreqRxCorrected2(34:59) ./ lambdas(27:end);

    
%     fineWindow(d)
%     doInnerProduct(ltsFreqRxCorrected2,ltsFreqRxCorrected2)
    innerProducts(d) = doInnerProduct(ltsFreqRxCorrected2,ltsFreqRxCorrected2);
    
    if fineWindow(d) == 1001
        ltsFreqRxCorrected2;
    end
    
    totalErrors(d) = sum(abs(ltsFreq-ltsFreqRxCorrected2));
%     totalErrors(d) = doInnerProduct([ltsFreqRxCorrected1; ltsFreqRxCorrected2], [ltsFreq; ltsFreq]);
end

quantileError = quantile(totalErrors, [0.3, 1])
% filteredErrors = totalErrors(totalErrors <= medianError)
fineTauTest = 0;
for d = 1:length(totalErrors)
   if totalErrors(d) <= quantileError
       fineTauTest = fineWindow(d);
   end
end

figure;
plot(fineWindow, abs(totalErrors),'.')
% hold on;
% plot(fineWindow, filteredErrors, 'r.')
title('Total error on LTS vs offset');

% errorDifferenceWithPrevious = 
% 
% for d = 0:length(totalErrors)-1;
%    
%     if totalErrors()
% 
% end
    
[~, pos] = min(totalErrors);
% [~, pos] = max(totalErrors);
fineTauMin = fineWindow(pos)
fineTau = fineTauTest


% Fine timing estimation
% coarseWindowLength = 60; % length of ISI 
% frameOffsetStartTwoLts = 10*length(stsTime) + twoLtsCpLength;
% fineWindow = (coarseTau-coarseWindowLength:coarseTau);
% innerProducts = zeros(1,length(fineWindow));
% twoLtsTime = repmat(ltsTime, 2,1);
% normalization = sum(abs(twoLtsTime).^2);
% for d = 1:length(fineWindow)
%     twoLtsTimeRx = signal(fineWindow(d)+frameOffsetStartTwoLts:fineWindow(d)+frameOffsetStartTwoLts+length(twoLtsTime)-1);
%     innerProducts(d) = doInnerProduct(twoLtsTimeRx, twoLtsTime)/normalization;
% end
% 
% figure;
% subplot(2,1,1),plot(fineWindow, abs(innerProducts),'.')
% title('Inner-products with long training sequence (Mag)');
% subplot(2,1,2),plot(fineWindow, angle(innerProducts),'.')
% title('Inner-products with long training sequence (Phase)');

% Fine timing estimation
% coarseWindowLength = 145; 
% if decoder.upsample
%     coarseWindowLength = coarseWindowLength*decoder.USF;
% end
% frameOffsetStartTwoLts = 10*length(stsTime) + twoLtsCpLength;
% fineWindow = (coarseTau-coarseWindowLength+1:coarseTau+coarseWindowLength-1);
% innerProducts = zeros(1,length(fineWindow));
% twoLtsTime = repmat(ltsTime, 2,1);
% normalization = sum(abs(twoLtsTime).^2);
% for d = 1:length(fineWindow)
%     twoLtsTimeRx = signal(fineWindow(d)+frameOffsetStartTwoLts:fineWindow(d)+frameOffsetStartTwoLts+length(twoLtsTime)-1);
%     innerProducts(d) = doInnerProduct(twoLtsTimeRx, twoLtsTime)/normalization;
% end
% % 
% [~, pos] = max(abs(innerProducts));
% fineTau = fineWindow(pos)
% 
% figure;
% subplot(2,1,1),plot(fineWindow, abs(innerProducts),'.')
% title('Inner-products with long training sequence (Mag)');
% subplot(2,1,2),plot(fineWindow, angle(innerProducts),'.')
% title('Inner-products with long training sequence (Phase)');

% tau = coarseTau + decoder.timingOffset;
tau = fineTau + decoder.timingOffset

% tau = 2024

end

