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

% mode = 'twoBoardsRxTx';
mode = 'simulation';

% Decoder parameters
decoderConfig()

config.Fs = 0.5e6; % Hz
config.M = 4;
config.map = qammap(config.M);
config.numTotalCarriers = 64;
marginSubCarriers = 0.17*config.numTotalCarriers;
% config.numZerosTop = 1;
% config.numZerosBottom = 1;
config.numZerosTop = ceil(marginSubCarriers/2);
config.numZerosBottom = floor(marginSubCarriers/2);
config.dcSubcarrier = 1; % 0 if used for data, 1 means that it isn't used to send data
config.numPilots = 2;
config.NFFT = config.numTotalCarriers;
config.NFFT_TS = 64;
config.CPLength = 16;
config.twoLtsCpLength = 32;
% config.twoLtsCpLength = 16;
config.numOFDMSymbolsPerFrame = 20; % Excluding preamble, must be >= 2
config.numUsedCarriers = config.numTotalCarriers - config.numZerosTop - config.numZerosBottom - config.dcSubcarrier;
config.numDataCarriers = config.numUsedCarriers - config.numPilots;

% Display config
config 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Training sequence for preamble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Short Training Sequence (STS), it is already extended to 64 and scaled as
% in 802.11a
load('stsFreq.mat');

% Long Training Sequence (LTS), it is already extended to 64 and scaled as
% in 802.11a
load('ltsFreq.mat');

load('CAcodes.mat');
ca = satCAcodes(1,:)/10;

fourStsTime = ifft(fftshift(stsFreq), config.NFFT_TS);
stsTime = sqrt(13/6)*fourStsTime(1:16); % cancel normalizations
ltsTime = ifft(fftshift(ltsFreq), config.NFFT_TS);

if strcmp(mode, 'simulation') || strcmp(mode, 'oneBoardTx') || strcmp(mode, 'loopback') || strcmp(mode, 'twoBoardsRxTx')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Information symbols Generation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Number of useful data symbol to send
    numDataSymbols = (config.numOFDMSymbolsPerFrame-1)*config.numDataCarriers
    dataSymbolsDec = randi([0,config.M-1], numDataSymbols, 1)
    % load('dataSymbolsDec.mat');

    % Get data symbols from constellation
    dataSymbols = modulator(dataSymbolsDec, config.map);

    % Pilot OFDM symbols (excluding pilot subcarriers) for channel estimation
    % Randomization is needed to avoid high PAPR
    % TODO: Use pseudorandomization to be able to use with USRPs on separate
    % hosts !
    % Note: To estimate the noise variance for MMSE equalization, it is needed
    % that re(X) = im(X) where X is one symbol of the pilot OFDM symbol.
    pilotOfdmSymbol = config.map(randi([0,1], config.numDataCarriers, 1)+2).'
    % load('pilotOfdmSymbol.mat');

    % Pilot subcarriers to correct for Sampling Frequency Offset (SFO) 
    % Note: the first symbol correspond to the pilot OFDM symbol
    pilotsSc1 = [config.map(1), repmat(config.map(config.M), 1, config.numOFDMSymbolsPerFrame-1)];
    pilotsSc2 = [config.map(4), repmat(config.map(config.M), 1, config.numOFDMSymbolsPerFrame-1)];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Transmitter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build the OFDM frame (preamble, pilots, data)
    [signalTx, dataFrame] = OFDMTransmitter(config, mode, dataSymbols, stsTime, ltsTime, pilotOfdmSymbol, pilotsSc1, pilotsSc2,ca);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;
% plot(abs(signalTx));
% xlabel('Samples');
% ylabel('Magnitude');
% title('Absolute value of transmitted signal');

% Transmission / Reception of the signal
if strcmp(mode, 'simulation') % Use simulator
    [signalRx, Fs] = simulatorOFDM(signalTx, config);
else % Use USRPs
    [signalRx,Fs] = rxTxUSRP(signalTx.', mode, config); 
    signalRx = signalRx.';
    signalRx = signalRx ./ max(abs(signalRx)); %TMP (?)
    
    % Truncate begining of the signal: contains garbage (how come ?)
    signalRx = signalRx(1e5:end);
end

if strcmp(mode, 'simulation') || strcmp(mode, 'twoBoardsRxTx') || strcmp(mode, 'oneBoardRx') || strcmp(mode, 'loopback')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Receiver
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmp(mode, 'oneBoardRx')
        dataFrame = [];
    end
    
    % Signal equalization and decoding
    dataSymbolsRx = OFDMReceiver(signalRx, config, stsTime, ltsTime, pilotOfdmSymbol, pilotsSc1, pilotsSc2, dataFrame, ca);
    
    % Demodulates symbols estimate(Hard Decision)
    dataSymbolsDecEst = demodulator(dataSymbolsRx, config.map);
    
    % Compute Symbol Error Rate (SER)
    SER = sum(dataSymbolsDec ~= dataSymbolsDecEst)/length(dataSymbolsDec)
end

fprintf('Instance terminated on %s \n\n',datestr(now));