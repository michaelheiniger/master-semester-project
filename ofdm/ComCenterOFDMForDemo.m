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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuration of data source

% bitsToSend = getRandomBits(10000);
bitsToSend = textFileToBits('textfile.txt');

bitsToSend = bitsToSend(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System configuration 
sc.Fs = 0.5e6; % [Hz], sampling frequency
sc.M = 4; % Size of M-QAM constellation
sc.map = qammap(sc.M); % M-QAM constellation vector
sc.numTotalCarriers = 64; % Total number of subcarriers, guard bands included (FFT size)
sc.numPilots = 2; % Number of subcarriers used as pilots (NOTE: system won't adapt to change)
sc.CPLength = 16; % Length of the cyclic prefix of the pilot and regular OFDM symbols (16 in IEEE 802.11a)
sc.twoLtsCpLength = 32; % Length of the cyclic prefix of the two LTSs of the preamble (32 in IEEE 802.11a)
sc.numOFDMSymbolsPerFrame = 100; % Excluding preamble, must be >= 2
sc.zeroFreqSubcarrier = 1; % 0 if used for data, 1 means that it is set to zero
marginSubCarriers = 0.17*sc.numTotalCarriers; % Ratio of subcarriers that are considered to be at the margin of the spectrum used
sc.numZerosTop = ceil(marginSubCarriers/2); % Amount of top (negative freq) subcarriers used as guard bands
sc.numZerosBottom = floor(marginSubCarriers/2); % Amount of bottom (positive freq) subcarriers used as guard bands
% config.numZerosTop = 1; > 0 since channel estimation is done on guard bands
% config.numZerosBottom = 1; > 0 since channel estimation is done on guard bands
sc.NFFT = sc.numTotalCarriers; % FFT size
sc.numUsedCarriers = sc.numTotalCarriers - sc.numZerosTop - sc.numZerosBottom - sc.zeroFreqSubcarrier; % number of subcarriers used to carry symbols (information or pilot)
sc.numDataCarriers = sc.numUsedCarriers - sc.numPilots; % number of subcarriers used to carry information symbols (i.e. not pilots)
sc.signalFieldLength = 64;
sc.signalFieldCpLength = 16;
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

% Configuration of OFDM receiver 1
rc1.timingAndFrequencyOffsetMethod = 'stsLtsOfdmDemod';
rc1.cfoCorrection = 1; % 1 if CFO should be corrected
rc1.cfoTracking = 1; % 1 if residual CFO should be tracked
rc1.sfoCorrection = 0; % 1 if SFO should be corrected
rc1.equalization = 1; % 1 if channel equalization should be performed

rc1.timingOffset = 0; % add offset to timing estimate, in number of samples.
rc1.upsample = 0; % 1 if the received signal should be upsample for timing synchronization
rc1.USF = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Training sequence for preamble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% C/A code 
ca = getCA()/10; % power reduction

% Short Training Sequence of IEEE 802.11a
[stsTime, stsFreq] = getSTS();
% Cancel IEEE 802.11a normalization                                                    TODO: explain why ...
stsTime = sqrt(13/6)*stsTime;
stsFreq = sqrt(13/6)*stsFreq;

% Long Training Sequence of IEEE 802.11a
[ltsTime, ltsFreq] = getLTS();

% Pilot OFDM symbols (excluding pilot subcarriers) for channel estimation
% Randomization is needed to avoid high PAPR
pilotOfdmSymbol = buildPilotOfdmSymbol(sc.numUsedCarriers, sc.map, sc.M);

% Pilot subcarriers for CFO tracking and SFO correction
% Note: the first OFDM symbol correspond to the pilot OFDM symbol
pilotSubcarrier1 = repmat(pilotOfdmSymbol(1), 1, sc.numOFDMSymbolsPerFrame-1);
pilotSubcarrier2 = repmat(pilotOfdmSymbol(end), 1, sc.numOFDMSymbolsPerFrame-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmitter
% Fetch the bits to send and builds the OFDM frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isInstanceTransmitter(mode)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Information symbols generation
   
    numBitsToSend = length(bitsToSend);
    numInfoSymbolsNeeded = ceil(numBitsToSend / log2(sc.M));
            
    numOFDMSymbolsNeeded = ceil(numInfoSymbolsNeeded/sc.numDataCarriers);
    if numOFDMSymbolsNeeded > sc.numOFDMSymbolsPerFrame-1
        error('Number of OFDM symbols per frame exceeded')
    end
    
    numBitsPadding = (sc.numOFDMSymbolsPerFrame-1)*sc.numDataCarriers*log2(sc.M)-numBitsToSend;
        
    bitsToSendWithPadding = [bitsToSend; randi([0,1], numBitsPadding, 1)];
    
    % Transform bits into information symbols from M-QAM constellation
    decInfoSymbols = bitsToDecSymbols(bitsToSendWithPadding, sc.M);
    infoSymbols = modulator(decInfoSymbols, qammap(sc.M));
    
    % Signal field uses BPSK (+/-1)
    signalSymbols = [de2bi(numBitsToSend, sc.numBitsForPayloadSize).';...
        randi([0,1], sc.numUsedCarriers-sc.numBitsForPayloadSize, 1)]*(-2)+1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OFDM transmitter
    % Build the OFDM frames and serialize them into signalTx
    % dataFrame is the frame containing the subcarriers used to carry
    % symbols (information or pilot) and is returned only for comparison
    % with received version (see OFDMReceiver)
    [signalTx, dataFrame] = OFDMTransmitter(sc, infoSymbols, stsTime, ltsTime, signalSymbols, pilotOfdmSymbol, pilotSubcarrier1, pilotSubcarrier2, ca);
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
    [infoSymbolsRx, numUsefulBitsRx] = OFDMReceiver(sc, rc1, signalRx, stsTime, ltsTime, pilotOfdmSymbol, pilotSubcarrier1, pilotSubcarrier2, dataFrame, ca);
    bitsRx = getBitsFromReceivedSymbols(infoSymbolsRx, sc.map, sc.M);
   numUsefulBitsRx
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