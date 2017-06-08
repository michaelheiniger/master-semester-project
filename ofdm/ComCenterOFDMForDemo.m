clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit and receive bits using OFDM modulation
% (See the bottom of this file for acronyms and abbreviations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Instance started on %s \n\n', datestr(now))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controls and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Possible modes: 
% - 'simulation'        can simulate different impairments: AWGN, CFO, SFO,
%                       ...
% - 'loopback':         one USRP board receives its own transmission 
%                       (--> Tx and Rx ports need to be linked by cable)
% - 'oneBoardTx':       one USRP board transmits per MATLAB instance
% - 'oneBoardRx':       one USRP board receives per MATLAB instance
% - 'twoBoardsRxTx':    two USRP boards on the same MATLAB instance: 
%                       one transmits, the other receives
% Note: USRP boards are always used in half-duplex mode except for the
% loopback mode
mode = 'simulation';
% mode = 'twoBoardsRxTx';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration of data source (source of bits to send)
% Two modes:
% - random bits: randomly generated bits, parameter: numRandomBits
% - textFile, bits from text file, parameter: textFilePath

% dataSource = 'random';
% numRandomBits = 10000;
dataSource = 'textFile';
textFilePath = 'textfile.txt';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System configuration 
sc.Fs = 0.5e6; % [Hz], sampling frequency
sc.M = 4; % Size of M-QAM constellation
sc.map = qammap(sc.M); % M-QAM constellation vector
sc.numTotalCarriers = 64; % Total number of subcarriers, guard bands included (FFT size)
sc.numPilots = 2; % Number of subcarriers used as pilots (NOTE: system won't adapt to change)
sc.CPLength = 16; % Length of the cyclic prefix of the pilot and regular OFDM symbols (16 in IEEE 802.11a)
sc.twoLtsCpLength = 32; % Length of the cyclic prefix of the two LTSs of the preamble (32 in IEEE 802.11a)
sc.numOFDMSymbolsPerFrame = 60; % Excluding preamble, must be >= 3 (pilot OFDM symbol, signal OFDM symbol)
sc.zeroFreqSubcarrier = 1; % 0 if used for data, 1 means that it is set to zero
marginSubCarriers = 0.17*sc.numTotalCarriers; % Ratio of subcarriers that are considered to be at the margin of the spectrum used
% sc.numZerosTop = ceil(marginSubCarriers/2); % Amount of top (negative freq) subcarriers used as guard bands
% sc.numZerosBottom = floor(marginSubCarriers/2); % Amount of bottom (positive freq) subcarriers used as guard bands
sc.numZerosTop = 1; %> 0 since channel estimation is done on guard bands
sc.numZerosBottom = 1; %> 0 since channel estimation is done on guard bands
sc.NFFT = sc.numTotalCarriers; % FFT size
sc.numUsedCarriers = sc.numTotalCarriers - sc.numZerosTop - sc.numZerosBottom - sc.zeroFreqSubcarrier; % number of subcarriers used to carry symbols (information or pilot)
sc.numDataCarriers = sc.numUsedCarriers - sc.numPilots; % number of subcarriers used to carry information symbols (i.e. not pilots)
sc.numBitsForPayloadSize = 16; % first bits of signal field
sc % Display config

% The noise variance is estimated using the guard bands so there must be at
% least one.
if sc.numZerosBottom+sc.numZerosTop == 0
   error('There is no zero subcarriers: cannot measure the noise variance using MMSE !'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Receivers configuration
% Several receivers can be defined for comparison
% Methods for timing and frequency offset estimation are:
% - stsLtsOfdmDemod: use IEEE 802.11a preamble and demodulates preamble
%   symbols to compute estimates
% - caTimeDomain: use Gold sequence C/A code as a preamble and time-domain
%   channel estimation to compute estimates

% Configuration of OFDM receiver
rc.timingAndFrequencyOffsetMethod = 'stsLtsOfdmDemod';
% rc.timingAndFrequencyOffsetMethod = 'ideal';
% rc.manualTiming = 1001; % Set the timing sample regardless of the actual estimate
% rc.manualCFO = 456; % Set the timing sample regardless of the actual estimate
rc.cfoCorrection = 1; % 1 if CFO should be corrected
rc.cfoTracking = 1; % 1 if residual CFO should be tracked
rc.sfoCorrection = 0; % 1 if SFO should be corrected
rc.equalization = 1; % 1 if channel equalization should be performed

rc.timingOffset = 0; % add offset to timing estimate, in number of samples.
rc.upsample = 0; % 1 if the received signal should be upsample for timing synchronization
rc.USF = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmitter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isInstanceTransmitter(mode)
    
    % Fetch the bits to send
    if strcmp(dataSource, 'random')
        bitsToSend = getRandomBits(numRandomBits);    
    elseif strcmp(dataSource, 'textFile')
        bitsToSend = textFileToBits(textFilePath);
    else
        error('Unknown data source');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OFDM transmitter
    % Build the OFDM frame and serialize it into signalTx
    % dataFrame is the frame containing the subcarriers used to carry
    % symbols (information or pilot) and is returned only for comparison
    % with received version (see OFDMReceiver)
    [signalTx, dataFrame, bitsToSendWithPadding, infoSymbols] = OFDMTransmitter(bitsToSend, sc);
else
    % If the current MATLAB instance if NOT a transmitter, the signal to
    % transmit is empty (necessary for USRP code)
    signalTx = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmission / Reception of the signal over the channel
% Simulator or USRP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotSignalMagnitude(signalTx, 'Samples', 'Absolute value of transmitted signal')

% Transmission / Reception of the signal
if strcmp(mode, 'simulation') % Use simulator
    
    % Add garbage before and after OFDM frame and repeat
    signalTx = repmat([zeros(1000,1); signalTx; zeros(1000,1)],1,1);
    
    % Simulate the impairements of the channel on the signal
    signalRx = simulatorOFDM(signalTx, sc.Fs);
else % Use USRPs
    
    % Repeat the signal
    signalTx = repmat(signalTx,1000,1);
        
    % Transmit / Receive with USRPs
    signalRx = rxTxUSRP(signalTx, mode, sc.Fs);
    
    % Truncate beginning of the signal: contains garbage (why ?)
    signalRx = signalRx(1e5:end);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OFDM Receivers
% Feed the different OFDM receivers with the received signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isInstanceReceiver(mode)
    if strcmp(mode, 'oneBoardRx')
        dataFrame = [];
    end
    
    % Receiver
    [bitsRx, infoSymbolsRx, numUsefulBitsRx, ~] = OFDMReceiver(sc, rc, signalRx, dataFrame);
    bitsRx = getBitsFromReceivedSymbols(infoSymbolsRx, sc.map, sc.M);
    disp(['Number of useful bits sent (extracted from SIGNAL field): ' num2str(numUsefulBitsRx)]);
    
    % Remove paddind bits
    usefulBits = bitsRx(1:numUsefulBitsRx);
    
    if not(strcmp(mode, 'oneBoardRx'))
        % Compute Bit and Symbol Error Rates of all receivers
        [BER, SER] = symbolsAndBitsStats(infoSymbols, infoSymbolsRx, bitsToSendWithPadding, bitsRx, sc.map)
    end
    
    disp('Text decoded by receiver 1');
    textRx = bitsToText(usefulBits)
end

fprintf('Instance terminated on %s \n\n',datestr(now));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acronyms and abbreviations:
% AWGN: Additive White Gaussian Noise
% CFO: Carrier Frequency Offset
% SFO: Sampling Frequency Offset
% USRP: Universal Software Radio Peripheral
% Tx: Transmitter
% Rx: Receiver
% CP: Cyclic Prefix
% STS: Short Training Sequence
% LTS: Long Training Sequence
% CA: Coarse Acquisition codes (in reference to GPS C/A codes)
% USF: Up-Sampling Factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%