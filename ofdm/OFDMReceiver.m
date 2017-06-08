function [bitsRx, infoSymbolsRx, numUsefulBitsRx, timingEst] = OFDMReceiver(systemConfig, receiverConfig, signalRx, dataFrame)
%OFDMRECEIVER Summary of this function goes here
% Receiver is configured by struct receiverConfig
% Parameters of the receiver are:
% - timingOffset = 0; % add offset to timing estimate, in number of samples.
% - timingOffsetCorrection = 0; % 1 if timing offset should be corrected
% - cfoCorrection = 1; % 1 if CFO should be corrected
% - cfoTracking = 1; % 1 if residual CFO should be tracked
% - sfoCorrection = 0; % 1 if SFO should be corrected
% - equalization = 1; % 1 if channel equalization should be performed
% - upsample = 0; % 1 if the received signal should be upsample for timing synchronization
% - USF = 10;

sc = systemConfig;
rc = receiverConfig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fetch training sequences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get C/A code 
ca = getCA()/10; % power reduction

% Get Short Training Sequence of IEEE 802.11a
[stsTime, ~] = getSTS();

% Get Long Training Sequence of IEEE 802.11a
[ltsTime, ~] = getLTS();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fetch OFDM pilot symbols and pilot subcarriers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pilot OFDM symbols (excluding pilot subcarriers) for channel estimation
% Randomization is needed to avoid high PAPR
pilotOfdmSymbol = buildPilotOfdmSymbol(sc.numUsedCarriers, sc.map, sc.M);

% Pilot subcarriers for CFO tracking and SFO correction
% Note: the first OFDM symbol correspond to the pilot OFDM symbol
pilotSubcarrier1 = repmat(pilotOfdmSymbol(1), 1, sc.numOFDMSymbolsPerFrame-1);
pilotSubcarrier2 = repmat(pilotOfdmSymbol(end), 1, sc.numOFDMSymbolsPerFrame-1);


plotSignalMagnitude(signalRx, 'Samples', 'Absolute value of received signal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OFDM demodulation of the received signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Timing synchronization and Frequency offset estimation and correction
[frameRx, timingEst] = timingAndFrequencyOffsetCorrection(systemConfig, receiverConfig, signalRx, stsTime, ltsTime, ca);

% Remove preamble (1 CA, 10 STS, CP, 2 LTS, CP, Signal Field)
preambleLength = length(ca)+10*length(stsTime)+sc.twoLtsCpLength+2*length(ltsTime);
dataRxCorrected = frameRx(1+preambleLength:end);

% Reshape into a matrix and remove cyclic prefix
dataRxIfftWithCp = reshape(dataRxCorrected, sc.CPLength+sc.numTotalCarriers, sc.numOFDMSymbolsPerFrame);
dataRxIFFT = dataRxIfftWithCp(sc.CPLength+1:end,:);

% FFT: ifftshift() is needed to get the subcarriers in increasing order of
% indices (i.e. negative SCs at the top, positive SCs at the bottom of the
% frame)
ofdmSymbolsRx = ifftshift(fft(dataRxIFFT, sc.NFFT),1);

% Extract received pilot OFDM symbol
pilotOfdmSymbolRx = ofdmSymbolsRx(:,1);

% Measure noise variance using guard bands (zero subcarriers)
zeroSubcarriersTop = pilotOfdmSymbolRx(1:sc.numZerosTop);
zeroSubcarriersBottom = pilotOfdmSymbolRx(end-sc.numZerosBottom+1);
noiseVariance = 1/sc.numTotalCarriers*var([zeroSubcarriersTop; zeroSubcarriersBottom]); % Normalization since noise is added in time domain

% Remove guard bands (i.e. outer-zero-subcarriers)
ofdmSymbolsRx = ofdmSymbolsRx(1+sc.numZerosTop:end-sc.numZerosBottom,:);
pilotOfdmSymbolRx = pilotOfdmSymbolRx(1+sc.numZerosTop:end-sc.numZerosBottom);

% Remove DC subcarrier (i.e. k = 0 in FFT) if present 
if sc.zeroFreqSubcarrier
    ofdmSymbolsRx = [ofdmSymbolsRx(1:sc.numTotalCarriers/2-sc.numZerosTop,:);...
                     ofdmSymbolsRx(2+sc.numTotalCarriers/2-sc.numZerosTop:end,:)];
    
    pilotOfdmSymbolRx = [pilotOfdmSymbolRx(1:sc.numTotalCarriers/2-sc.numZerosTop);...
                         pilotOfdmSymbolRx(2+sc.numTotalCarriers/2-sc.numZerosTop:end)];
end

% Channel estimation using Minimum Mean Squared Error (MMSE)
lambdasMMSE = channelEstimation(ofdmSymbolsRx(:,2:end), pilotOfdmSymbolRx, pilotOfdmSymbol, sc.numUsedCarriers, noiseVariance);

% Extract the most faded subcarriers indices
lowestChannelsCoeffsIndices = extractMostFadedSubcarriersIndices(lambdasMMSE, 0.1, sc.numUsedCarriers);

% Simple channel estimation
lambdas = pilotOfdmSymbolRx ./ pilotOfdmSymbol;

% Get indices of used subcarriers (e.g -26,-25,...,+26)
if sc.zeroFreqSubcarrier
    usedSubcarriersIndices = [-sc.numTotalCarriers/2+sc.numZerosTop:-1, 1:sc.numTotalCarriers/2-1-sc.numZerosBottom];
else
    usedSubcarriersIndices = -sc.numTotalCarriers/2+sc.numZerosTop:sc.numTotalCarriers/2-1-sc.numZerosBottom;
end

% Plot of the channel coefficients
figure;
xAxis = usedSubcarriersIndices;

subplot(1,2,1),
plot(xAxis, abs(lambdas),'r.-');
hold on;
plot(xAxis, abs(lambdasMMSE), 'b.-');
ylabel('Magnitude');
xlabel('Subcarriers index');
legend('Division','MMSE');
grid on;

subplot(1,2,2),
plot(xAxis, angle(lambdas),'r.-');
hold on;
plot(xAxis, angle(lambdasMMSE), 'b.-');
ylabel('Phase');
xlabel('Subcarriers index');
legend('Division','MMSE');
grid on;
axes;
set(gca,'Visible','off');
set(title('Channel coefficient estimates'),'Visible','on'); 

% Save received symbols to compare with symbols corrected for SFO and
% equalized
ofdmSymbolsRxCorrected = ofdmSymbolsRx;

plotChannelCoeffFrame(systemConfig, ofdmSymbolsRxCorrected, dataFrame, 'Before equalization');

% Equalization of the symbols
if rc.equalization
    ofdmSymbolsRxCorrected = ofdmSymbolsRxCorrected ./ repmat(lambdasMMSE, 1, sc.numOFDMSymbolsPerFrame);
%     ofdmSymbolsRxCorrected = ofdmS ymbolsRxCorrected ./ repmat(lambdas, 1, c.numOFDMSymbolsPerFrame);

    plotChannelCoeffFrame(systemConfig, ofdmSymbolsRxCorrected, dataFrame, 'After equalization, before CFO tracking');
end

% Complete the pilots with the OFDM pilot symbol for CFO trakcing and SFO
% correction
pilotSubcarrier1 = [pilotOfdmSymbol(1), pilotSubcarrier1];
pilotSubcarrier2 = [pilotOfdmSymbol(end), pilotSubcarrier2];

if rc.cfoTracking
    % Track and remove the residual drifting CFO present in all OFDM symbols due to
    % imperfect estimate of the CFO
    ofdmSymbolsRxCorrected = carrierFrequencyOffsetTracking(ofdmSymbolsRxCorrected, pilotSubcarrier1, pilotSubcarrier2);
    
    plotChannelCoeffFrame(systemConfig, ofdmSymbolsRxCorrected, dataFrame, 'After CFO tracking, before SFO');
end

if rc.sfoCorrection
    % Sampling Frequency Offset correction
    ofdmSymbolsRxCorrected = samplingFrequencyOffsetCorrection(ofdmSymbolsRxCorrected, sc, pilotSubcarrier1, pilotSubcarrier2);
    
    plotChannelCoeffFrame(systemConfig, ofdmSymbolsRxCorrected, dataFrame, 'After SFO');
end

% Extract SIGNAL OFDM symbol (without pilots)
signalFieldOfdmSymbolRx = ofdmSymbolsRxCorrected(2:end-1,2);

figure;
plot(signalFieldOfdmSymbolRx,'.');
hold on;
plot([-1 1], [0 0], 'gx');
xlabel('In-phase');
ylabel('Quadrature');
title('BPSK SIGNAL symbols')

% Remove SIGNAL OFDM symbol from receiver symbols
ofdmSymbolsRxCorrected = ofdmSymbolsRxCorrected(:,[1,3:end]);

% Signal field bits extraction
signalFieldOfdmSymbolRx(real(signalFieldOfdmSymbolRx) >= 0) = 0;
signalFieldOfdmSymbolRx(real(signalFieldOfdmSymbolRx) < 0) = 1;
numUsefulBitsRxBin = signalFieldOfdmSymbolRx(1:sc.numBitsForPayloadSize);
numUsefulBitsRx = bi2de(numUsefulBitsRxBin(:).');

% Split symbols sent on outer-subcarriers and inner-subcarriers to show the
% effect of SFO on outer-subcarriers. Enabled only when no guard bands are
% used. Symbols sent on outer-subcarriers are plotted in red, others in
% blue
outerSubcarrierRatio = 0.19;
if sc.numZerosTop+sc.numZerosBottom == 0
    numOuterSubcarriers = ceil(sc.numTotalCarriers*outerSubcarrierRatio);
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
% figure;
% plot(infoSymbolOuterSCEst, 'r.');
% hold on;
% grid on;
% plot(infoSymbolInnerSCEst, 'b.');
% plot(qammap(c.M), 'gx');
% xlabel('In-phase');
% ylabel('Quadrature');
% title([num2str(c.M) '-QAM constellation at receiver']);
% legend(['Outer SC symbols (2*' num2str(outerSubcarrierRatio*100/2) '%)'], ['Inner SC symbols (' num2str(100*(1-outerSubcarrierRatio)) '%)'], 'QAM points');
% axis([-1,1,-1,1])
% axis([-1.5,1.5,-1.5,1.5]);

% Plot symbols on M-QAM constellation: symbols are corrected for CFO, SFO 
% and equalized. Pilot OFDM symbol and symbols sent over pilot subcarriers 
% are included.

figure;
plot(infoSymbolOuterSCCorrectedEst, 'r.');
hold on;
grid on;
plot(infoSymbolInnerSCCorrectedEst, 'b.');
plot(ofdmSymbolsRxCorrected(lowestChannelsCoeffsIndices,:),'m.');
plot(qammap(sc.M), 'gx');
xlabel('In-phase');
ylabel('Quadrature');
title([num2str(sc.M) '-QAM constellation at receiver']);
axis([-1.5,1.5,-1.5,1.5]);

% legend(['Outer SC symbols (2*' num2str(outerSubcarrierRatio*100/2) '%)'], ['Inner SC symbols (' num2str(100*(1-outerSubcarrierRatio)) '%)'], 'QAM points');

% Remove pilots (i.e. remove first and last subcarriers) 
ofdmSymbolsRxCorrected = ofdmSymbolsRxCorrected(2:end-1,:);

% Remove OFDM pilot symbol
ofdmSymbolsRxCorrected = ofdmSymbolsRxCorrected(:,2:end);

% Return vector of data symbols: demodulation is needed !
infoSymbolsRx = ofdmSymbolsRxCorrected(:);

% Demodulates information symbols into bits
bitsRx = getBitsFromReceivedSymbols(infoSymbolsRx, sc.map, sc.M);

end

