function [dataSymbolsRx] = OFDMReceiver(signalRx, config, stsTime, ltsTime, pilotOfdmSymbol, pilots1, pilots2, dataFrame)
%OFDMRECEIVER Summary of this function goes here
%   Detailed explanation goes here

global decoder;

c = config;

figure;
plot(abs(signalRx));
xlabel('Samples');
ylabel('Magnitude');
title('Absolute value of received signal');

% Timing synchronization and Frequency offset estimation
[frameRxCorrected, ~] = timingAndFreqOffsetCorrection(signalRx, stsTime, ltsTime, c.numTotalCarriers, c.CPLength, c.twoLtsCpLength, c.numOFDMSymbolsPerFrame, c.Fs);

% Remove preamble
preambleLength = 10*length(stsTime)+c.twoLtsCpLength+2*length(ltsTime);
dataRxCorrected = frameRxCorrected(1+preambleLength:end);

% Reshape into a matrix and remove cyclic prefix
dataRxIfftCP = reshape(dataRxCorrected, c.CPLength+c.numTotalCarriers, c.numOFDMSymbolsPerFrame);
dataRxIFFT = dataRxIfftCP(c.CPLength+1:end,:);

% FFT: ifftshift() is needed to get the subcarriers in increasing order of
% indices (i.e. negative SCs at the top, positive SCs at the bottom of the
% frame)
ofdmSymbolsRx = ifftshift(fft(dataRxIFFT, c.NFFT),1);

% Remove guard bands (i.e. outer-zero-subcarriers)
ofdmSymbolsRx = ofdmSymbolsRx(1+c.numZerosTop:end-c.numZerosBottom,:);

% Remove DC subcarrier (i.e. k = 0 in FFT) if present 
if c.dcSubcarrier
    ofdmSymbolsRx = [ofdmSymbolsRx(1:c.numTotalCarriers/2-c.numZerosTop,:);...
                     ofdmSymbolsRx(2+c.numTotalCarriers/2-c.numZerosTop:end,:)];
end

% Save received symbols to compare with symbols corrected for SFO and
% equalized
ofdmSymbolsRxCorrected = ofdmSymbolsRx;

% Add pilot subcarrier to OFDM pilot symbol
pilotOfdmSymbol = [pilots1(1); pilotOfdmSymbol; pilots2(1)];

% Correct for timing offset (not actually used with USRP, only to makes
% plots to illustrate what happens when timing offset occurs)
if decoder.timingOffsetCorrection
    ofdmSymbolsRxCorrected = timingOffsetCorrection(ofdmSymbolsRxCorrected, pilotOfdmSymbol);
end

% if decoder.sfoCorrection
%     % Sampling Frequency Offset correction
%     ofdmSymbolsRxCorrected = samplingFrequencyOffsetCorrection(ofdmSymbolsRxCorrected, pilots1, pilots2, c.numZerosTop, c.numZerosBottom, c.numTotalCarriers, c.dcSubcarrier);
% end
   
% Get indices of used subcarriers (e.g -26,-25,...,+26)
if c.dcSubcarrier
    usedSubcarriersIndices = [-c.numTotalCarriers/2+c.numZerosTop:-1, 1:c.numTotalCarriers/2-1-c.numZerosBottom];
else
    usedSubcarriersIndices = -c.numTotalCarriers/2+c.numZerosTop:c.numTotalCarriers/2-1-c.numZerosBottom;
end

% Channel estimation using Minimum Mean Squared Error (MMSE)
lambdasMMSE = channelEstimation(ofdmSymbolsRxCorrected, pilotOfdmSymbol, c.numUsedCarriers);
lambdas = ofdmSymbolsRxCorrected(:,1) ./ pilotOfdmSymbol;

figure;
xAxis = usedSubcarriersIndices;
subplot(1,2,1),
plot(xAxis, abs(lambdas),'r.-');
hold on;
plot(xAxis, abs(lambdasMMSE), 'b.-');
ylabel('Magnitude');
xlabel('Subcarriers index');
grid on;

subplot(1,2,2),
plot(xAxis, angle(lambdas),'r.-');
hold on;
plot(xAxis, angle(lambdasMMSE), 'b.-');
ylabel('Phase');
xlabel('Subcarriers index');
grid on;

% Plot of the channel coefficients
% figure;
% xAxis = usedSubcarriersIndices;
% plot(xAxis, abs(lambdas),'.-');
% ylabel('Magnitude');
% xlabel('Subcarriers index');
% title('Channel coefficients estimates');
% grid on;

% Equalization of the symbols
if decoder.equalization
    ofdmSymbolsRxCorrected = ofdmSymbolsRxCorrected ./ repmat(lambdasMMSE, 1, c.numOFDMSymbolsPerFrame);
%     ofdmSymbolsRxCorrected = ofdmSymbolsRxCorrected ./ repmat(lambdas, 1, c.numOFDMSymbolsPerFrame);
end

% Plot of the channel coefficients for the whole frame (using division and
% NOT MMSE)
if c.dcSubcarrier
    % Compute channel coefficients, DC subcarrier excluded (it is only zeros)
    lambdaAllSymbols = ofdmSymbolsRxCorrected./[dataFrame(1+c.numZerosTop:c.numTotalCarriers/2,:);dataFrame(2+c.numTotalCarriers/2:end-c.numZerosBottom,:)];
else    
    % Compute channel coefficients, DC subcarrier included
    lambdaAllSymbols = ofdmSymbolsRxCorrected./dataFrame(1+c.numZerosTop:end-c.numZerosBottom,:);
end
figure;
xAxis = usedSubcarriersIndices;
subplot(2,1,1),plot(xAxis, angle(lambdaAllSymbols),'.-');
ylabel('Phase [rad]');
xlabel('Subcarriers index');
title('Phase evolution over subcarriers for the whole frame');
grid on;
% axis([-35,35,-0.5,0.5]);
subplot(2,1,2),plot(xAxis, abs(lambdaAllSymbols),'.-');
ylabel('Magnitude');
xlabel('Subcarriers index');
title('Magnitude evolution over subcarriers for the whole frame');
grid on;
% axis([-35,35,0.8,1.1]);

if decoder.sfoCorrection
    % Sampling Frequency Offset correction
    ofdmSymbolsRxCorrected = samplingFrequencyOffsetCorrection(ofdmSymbolsRxCorrected, pilots1, pilots2, c.numZerosTop, c.numZerosBottom, c.numTotalCarriers, c.dcSubcarrier);
end

% Split symbols sent on outer-subcarriers and inner-subcarriers to show the
% effect of SFO on outer-subcarriers. Enabled only when no guard bands are
% used. Symbols sent on outer-subcarriers are plotted in red, others in
% blue
outerSubcarrierRatio = 0.19;
if c.numZerosTop+c.numZerosBottom == 0
    numOuterSubcarriers = ceil(c.numTotalCarriers*outerSubcarrierRatio);
    if mod(numOuterSubcarriers, 2) ~= 0
        numOuterSubcarriers = numOuterSubcarriers + 1;
    end
    disp(['Number of outer-subcarriers: ' num2str(numOuterSubcarriers)]);
    
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

% Plot symbols on M-QAM constellation: symbols are corrected for CFO but
% not SFO and not equalized. Pilot OFDM symbol and symbols sent over pilot
% subcarriers are included.
figure;
plot(infoSymbolOuterSCEst, 'r.');
hold on;
grid on;
plot(infoSymbolInnerSCEst, 'b.');
plot(qammap(c.M), 'gx');
xlabel('In-phase');
ylabel('Quadrature');
title([num2str(c.M) '-QAM constellation at receiver']);
% legend(['Outer SC symbols (2*' num2str(outerSubcarrierRatio*100/2) '%)'], ['Inner SC symbols (' num2str(100*(1-outerSubcarrierRatio)) '%)'], 'QAM points');
% axis([-1,1,-1,1])

% Plot symbols on M-QAM constellation: symbols are corrected for CFO, SFO 
% and equalized. Pilot OFDM symbol and symbols sent over pilot subcarriers 
% are included.
figure;
plot(infoSymbolOuterSCCorrectedEst, 'r.');
hold on;
grid on;
plot(infoSymbolInnerSCCorrectedEst, 'b.');
plot(qammap(c.M), 'gx');
xlabel('In-phase');
ylabel('Quadrature');
% title([num2str(c.M) '-QAM constellation at receiver (SFO corrected and equalized)']);
title([num2str(c.M) '-QAM constellation at receiver']);

% legend(['Outer SC symbols (2*' num2str(outerSubcarrierRatio*100/2) '%)'], ['Inner SC symbols (' num2str(100*(1-outerSubcarrierRatio)) '%)'], 'QAM points');

% Remove pilots (i.e. remove first and last subcarriers) 
ofdmSymbolsRxCorrected = ofdmSymbolsRxCorrected(2:end-1,:);

% Remove OFDM pilot symbol
ofdmSymbolsRxCorrected = ofdmSymbolsRxCorrected(:,2:end);

% Return vector of data symbols: demodulation is needed !
dataSymbolsRx = ofdmSymbolsRxCorrected(:);

end

