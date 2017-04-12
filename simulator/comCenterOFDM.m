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

config.M = 4;
config.map = qammap(config.M);
config.numTotalCarriers = 128;
config.numZeros = 13;
config.numPilots = 2;
config.NFFT = config.numTotalCarriers;
config.NFFT_TS = 64;
config.CPLength = 16;
config.preamble2CPLength = 32;
config.numOFDMSymbolsPerFrame = 30; % Excluding preamble
config.numDataCarriers = config.numTotalCarriers - 2*config.numZeros;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preambles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('shortTrainingSeq.mat');
load('longTrainingSeq.mat');

shortTrainingSeq = [zeros(6,1); sqrt(13/6)*shortTrainingSeq; zeros(5,1)];
shortTrainingSeqIFFT = ifft(fftshift(shortTrainingSeq), config.NFFT_TS/4);

longTrainingSeq = [zeros(6,1); longTrainingSeq; zeros(5,1)];
longTrainingSeqIFFT = ifft(fftshift(longTrainingSeq), config.NFFT_TS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information symbols Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numSymbols = (config.numOFDMSymbolsPerFrame-1)*(config.numDataCarriers-config.numPilots); % -1 to account for the pilot symbol
symbols = modulator(randi([0,config.M-1], numSymbols, 1), config.map);

pilotOfdmSymbol = modulator(randi([0,config.M-1], config.numDataCarriers, 1), config.map);

pilots1 = ones(1, config.numOFDMSymbolsPerFrame-1);
pilots2 = ones(1, config.numOFDMSymbolsPerFrame-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmitter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[signalTx, dataFrame] = OFDMTransmitter(config, mode, symbols, shortTrainingSeqIFFT, longTrainingSeqIFFT, pilotOfdmSymbol, pilots1, pilots2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(abs(signalTx));
title('Absolute value of transmitted signal');

if strcmp(mode, 'simulation')
    [signalRx, Fs] = simulatorOFDM(signalTx);
else
    signalTx = padWithZeros(signalTx, samplesPerFrame);
    [signalRx,Fs] = rxTxUSRP(signalTx.', samplesPerFrame, mode); 
    signalRx = signalRx.'; 
    signalRx = signalRx(1e5:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(mode, 'simulation') || strcmp(mode, 'twoBoardsRxTx') || strcmp(mode, 'oneBoardRx') || strcmp(mode, 'loopback')
    symbolsRx = OFDMReceiver(signalRx, config, shortTrainingSeqIFFT, longTrainingSeqIFFT, pilotOfdmSymbol, pilots1, pilots2, dataFrame);
end

fprintf('Instance terminated on %s \n\n',datestr(now));