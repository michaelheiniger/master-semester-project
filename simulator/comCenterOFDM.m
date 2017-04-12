clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to send and receive data through simulator or USRP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Instance started on %s \n\n',datestr(now))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controls and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Possible modes: 
% - 'simulation' The downside of DA pre-FFT synchronization is reduced transmission efficien:       simulates impairements such as AWGN, Doppler, clock offset,... 
% - 'loopback':         one USRP board receives its own transmission 
%                       --> Tx and Rx ports need to be linked by cable 
% - 'oneBoardTx':       one USRP board transmits 
% - 'oneBoardRx':       one USRP board receives
% - 'twoBoardsRxTx':    two USRP boards on the same computer: one transmits, the
%                       other receives
% USRP boards are always used in half-duplex mode

mode = 'simulation';
samplesPerFrame = 5e4;

% Decoder parameters
decoderConfig()

M = 4;
map = qammap(M);
numTotalCarriers = 128;
numZeros = 13;
numPilots = 2;
NFFT = numTotalCarriers;
NFFT_TS = 64;
CPLength = 16;
preamble2CPLength = 32;
numOFDMSymbolsPerFrame = 10; % Excluding preamble
numDataCarriers = numTotalCarriers - 2*numZeros;

% pilotsPosition = [1, numDataCarriers];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmitter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Information symbols generation
numSymbols = (numOFDMSymbolsPerFrame-1)*(numDataCarriers-numPilots); % -1 to account for the pilot symbol

symbols = modulator(randi([0,M-1], numSymbols, 1), map);
pilotOfdmSymbol = modulator(randi([0,M-1], numDataCarriers, 1), map);

% First preamble construction
load('shortTrainingSeq.mat');
shortTrainingSeq = [zeros(6,1); sqrt(13/6)*shortTrainingSeq; zeros(5,1)];
shortTrainingSeqIFFT = ifft(fftshift(shortTrainingSeq), NFFT_TS/4);
preamble1 = repmat(shortTrainingSeqIFFT, 10, 1);

% Second preamble construction
load('longTrainingSeq.mat');
longTrainingSeq = [zeros(6,1); longTrainingSeq; zeros(5,1)];
longTrainingSeqIFFT = ifft(fftshift(longTrainingSeq), NFFT_TS);
preamble2 = repmat(longTrainingSeqIFFT, 2, 1);

% Add CP to second preamble
secondPreambleIFFTCP = [preamble2(end-preamble2CPLength+1:end); preamble2];

% OFDM data frame construction
dataFrame = reshape(symbols, numDataCarriers-numPilots, numOFDMSymbolsPerFrame-1);

% Add pilot subcarriers (for post-FFT frequency offset correction due to
% sampling frequency offset)
pilots1 = ones(1, numOFDMSymbolsPerFrame-1);
pilots2 = ones(1, numOFDMSymbolsPerFrame-1);
dataFrame = [pilots1;...
            dataFrame; ...
            pilots2];

% Add pilot ofdm symbol (for channel estimation)
dataFrame = [pilotOfdmSymbol, dataFrame];

dataFrame = [zeros(numZeros, numOFDMSymbolsPerFrame); ...
            dataFrame; ...
            zeros(numZeros, numOFDMSymbolsPerFrame)];
            
% Apply ffshift to account for the negative frequencies since the MATLAB
% function ifft takes k=0,...,N-1 and not k=-N/2,...,0,...,N/2-1
% Apply ifft to get the time domain signal
dataFrameIFFT = ifft(fftshift(dataFrame,1), NFFT);
% Add CP to data frame
dataFrameIFFTCP = [dataFrameIFFT(end-CPLength+1:end,:);...
                   dataFrameIFFT];

% Serialization
signalTx = [preamble1;...
          secondPreambleIFFTCP;...
          dataFrameIFFTCP(:)];

if strcmp(mode, 'simulation')
    signalTx = repmat([zeros(1000,1); signalTx; zeros(1000,1)],10,1);
else 
    signalTx = repmat([zeros(1000,1); signalTx; zeros(1000,1)],100,1);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(abs(signalTx));
title('Absolute value of transmitted signal');

if strcmp(mode, 'simulation')
    [signalRx, Fs] = simulatorOFDM(signalTx);
else
    signalTx = padWithZeros(signalTx,samplesPerFrame);
    [signalRx,Fs] = rxTxUSRP(signalTx.', samplesPerFrame, mode); 
    signalRx = signalRx.'; 
    signalRx = signalRx(1e5:end);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode, 'simulation') || strcmp(mode, 'twoBoardsRxTx') || strcmp(mode, 'oneBoardRx') || strcmp(mode, 'loopback')
    figure;
    hax2=axes;
    plot(abs(signalRx));
    title('Absolute value of received signal');
    
    STSIFFTLength = length(shortTrainingSeqIFFT);
    LTSIFFTLength = length(longTrainingSeqIFFT);
    
    preamble1Length = 10*STSIFFTLength;
    preamble2Length = 2*LTSIFFTLength;
    
    % Timing synchronization and Frequency offset estimation
    [frameRxCorrected, frequencyOffsetEst] = timingAndFreqOffsetCorrection(signalRx, shortTrainingSeqIFFT, longTrainingSeqIFFT, numTotalCarriers, CPLength, preamble2CPLength,numOFDMSymbolsPerFrame);
    
    % Remove preamble1 and preamble2
    dataRxCorrected = frameRxCorrected(1+preamble1Length+preamble2CPLength+preamble2Length:end);
    
    % Reshape into a matrix and remove CP
    dataRxIfftCP = reshape(dataRxCorrected, CPLength+numTotalCarriers, numOFDMSymbolsPerFrame);
    dataRxIFFT = dataRxIfftCP(CPLength+1:end,:);
    
    % FFT
    ofdmSymbolsRx = ifftshift(fft(dataRxIFFT, NFFT),1);
    
    % Remove guard bands (i.e. zero-subcarriers)
    ofdmSymbolsRx = ofdmSymbolsRx(1+numZeros:end-numZeros,:);
        
    % Channel estimation
    channelCoefficients = channelEstimation(ofdmSymbolsRx, pilotOfdmSymbol, numTotalCarriers, numZeros);
    
    % Equalization
    ofdmSymbolsRxCorrected = ofdmSymbolsRx ./ repmat(channelCoefficients, 1, numOFDMSymbolsPerFrame);
    
    % Post FFT frequency offset estimaton and correction
    ofdmSymbolsRxCorrected = postFFTFrequencyOffsetCorrection(ofdmSymbolsRxCorrected, channelCoefficients, pilots1, pilots2, numDataCarriers, numZeros);
    
    % TEST
    figure;
    channelCoefficientsWholeFrame = ofdmSymbolsRx./dataFrame(1+numZeros:end-numZeros,:);
    subplot(2,1,1),plot(angle(channelCoefficientsWholeFrame),'.-');
    title('Phase of channel coefficients');
    subplot(2,1,2),plot(abs(channelCoefficientsWholeFrame),'.-');
    title('Magnitude of channel coefficients');
    grid on;

    % /TEST      
    
    outerSubcarrierRatio = 0.19;
    if numZeros == 0   
        numOuterSubcarriers = ceil(numTotalCarriers*outerSubcarrierRatio);
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
    plot(qammap(M), 'gx');
    title([num2str(M) '-QAM constellation at receiver']);
    legend(['Outer SC symbols (2*' num2str(outerSubcarrierRatio*100/2) '%)'], ['Inner SC symbols (' num2str(100*(1-outerSubcarrierRatio)) '%)'], 'QAM points');
    
    figure;
    plot(infoSymbolOuterSCCorrectedEst, 'r.');
    hold on;
    grid on;
    plot(infoSymbolInnerSCCorrectedEst, 'b.');
    plot(qammap(M), 'gx');
    title([num2str(M) '-QAM constellation at receiver (Corrected)']);
    legend(['Outer SC symbols (2*' num2str(outerSubcarrierRatio*100/2) '%)'], ['Inner SC symbols (' num2str(100*(1-outerSubcarrierRatio)) '%)'], 'QAM points');
    
end

fprintf('Instance terminated on %s \n\n',datestr(now));