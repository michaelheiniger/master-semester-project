function [symbolsRx] = OFDMReceiver(signalRx, config, shortTrainingSeqIFFT, longTrainingSeqIFFT, pilotOfdmSymbol, pilots1, pilots2, dataFrame)
%OFDMRECEIVER Summary of this function goes here
%   Detailed explanation goes here

c = config;

figure;
plot(abs(signalRx));
title('Absolute value of received signal');

% Timing synchronization and Frequency offset estimation
[frameRxCorrected, ~] = timingAndFreqOffsetCorrection(signalRx, shortTrainingSeqIFFT, longTrainingSeqIFFT, c.numTotalCarriers, c.CPLength, c.preamble2CPLength, c.numOFDMSymbolsPerFrame);

% Remove preamble1 and preamble2
preamble1Length = 10*length(shortTrainingSeqIFFT);
preamble2Length = 2*length(longTrainingSeqIFFT);
dataRxCorrected = frameRxCorrected(1+preamble1Length+c.preamble2CPLength+preamble2Length:end);

% Reshape into a matrix and remove CP
dataRxIfftCP = reshape(dataRxCorrected, c.CPLength+c.numTotalCarriers, c.numOFDMSymbolsPerFrame);
dataRxIFFT = dataRxIfftCP(c.CPLength+1:end,:);

% FFT
ofdmSymbolsRx = ifftshift(fft(dataRxIFFT, c.NFFT),1);

% Remove guard bands (i.e. zero-subcarriers)
ofdmSymbolsRx = ofdmSymbolsRx(1+c.numZeros:end-c.numZeros,:);

% Channel estimation
channelCoefficients = channelEstimation(ofdmSymbolsRx, pilotOfdmSymbol, c.numTotalCarriers, c.numZeros);

% Equalization
ofdmSymbolsRxCorrected = ofdmSymbolsRx ./ repmat(channelCoefficients, 1, c.numOFDMSymbolsPerFrame);

% Post FFT frequency offset estimaton and correction
ofdmSymbolsRxCorrected = postFFTFrequencyOffsetCorrection(ofdmSymbolsRxCorrected, channelCoefficients, pilots1, pilots2, c.numDataCarriers, c.numZeros);

% TEST
figure;
channelCoefficientsWholeFrame = ofdmSymbolsRx./dataFrame(1+c.numZeros:end-c.numZeros,:);
subplot(2,1,1),plot(angle(channelCoefficientsWholeFrame),'.-');
title('Phase of channel coefficients');
subplot(2,1,2),plot(abs(channelCoefficientsWholeFrame),'.-');
title('Magnitude of channel coefficients');
grid on;

% /TEST

outerSubcarrierRatio = 0.19;
if c.numZeros == 0
    numOuterSubcarriers = ceil(c.numTotalCarriers*outerSubcarrierRatio);
    if mod(numOuterSubcarriers, 2) ~= 0
        numOuterSubcarriers = numOuterSubcarriers + 1;
    end
    disp(['num outer subcarriers: ' num2str(numOuterSubcarriers)]);
    
    infoSymbolOuterSCEst = [ofdmSymbolsRx(1:numOuterSubcarriers/2,:); ofdmSymbolsRx(end-numOuterSubcarriers/2+1:end,:)];
    infoSymbolInnerSCEst = ofdmSymbolsRx(1+numOuterSubcarriers/2:end-numOuterSubcarriers/2,:);
    
    infoSymbolOuterSCCorrectedEst = [ofdmSymbolsRxCorrected(1:numOuterSubcarriers/2,:); ofdmSymbolsRxCorrected(end-numOuterSubcarriers/2+1:end,:)];
    infoSymbolInnerSCCorrectedEst = ofdmSymbolsRxCorrected(1+numOuterSubcarriers/2:end-numOuterSubcarriers/2,:);
else
    infoSymbolOuterSCEst = [];
    infoSymbolInnerSCEst = ofdmSymbolsRx;
    infoSymbolOuterSCCorrectedEst = [];
    infoSymbolInnerSCCorrectedEst = ofdmSymbolsRxCorrected;
end

figure;
plot(infoSymbolOuterSCEst, 'r.');
hold on;
grid on;
plot(infoSymbolInnerSCEst, 'b.');
plot(qammap(c.M), 'gx');
title([num2str(c.M) '-QAM constellation at receiver']);
legend(['Outer SC symbols (2*' num2str(outerSubcarrierRatio*100/2) '%)'], ['Inner SC symbols (' num2str(100*(1-outerSubcarrierRatio)) '%)'], 'QAM points');

figure;
plot(infoSymbolOuterSCCorrectedEst, 'r.');
hold on;
grid on;
plot(infoSymbolInnerSCCorrectedEst, 'b.');
plot(qammap(c.M), 'gx');
title([num2str(c.M) '-QAM constellation at receiver (Corrected)']);
legend(['Outer SC symbols (2*' num2str(outerSubcarrierRatio*100/2) '%)'], ['Inner SC symbols (' num2str(100*(1-outerSubcarrierRatio)) '%)'], 'QAM points');

symbolsRx = infoSymbolOuterSCCorrectedEst;

end

