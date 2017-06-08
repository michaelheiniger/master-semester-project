% dataFrame is returned for experiment purposes only (plots, ...)
function [signalTx, dataFrame, bitsToSendWithPadding, infoSymbols] = OFDMTransmitter(bitsToSend, systemConfig)

sc = systemConfig;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preamble building
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get C/A code 
ca = getCA()/10; % power reduction

% Get Short Training Sequence of IEEE 802.11a
[stsTime, ~] = getSTS();
% Cancel IEEE 802.11a normalization                                                    TODO: explain why ...
stsTime = sqrt(13/6)*stsTime; 

% Get Long Training Sequence of IEEE 802.11a
[ltsTime, ~] = getLTS();

% Build 10 STSs vector
tenSts = repmat(stsTime, 10, 1);

% Build 2 LTSs vector
twoLts = repmat(ltsTime, 2, 1);

% Add cyclic prefix to the two LTS
twoLtsWithCp = [twoLts(end-sc.twoLtsCpLength+1:end); twoLts];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bits modulation into information symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numBitsToSend = length(bitsToSend);
numInfoSymbolsNeeded = ceil(numBitsToSend / log2(sc.M));
            
numOFDMSymbolsNeeded = ceil(numInfoSymbolsNeeded/sc.numDataCarriers);
if numOFDMSymbolsNeeded > sc.numOFDMSymbolsPerFrame-1
    error('Number of OFDM symbols per frame exceeded')
end

% "-2" because of pilot and signal OFDM symbols
numBitsPadding = (sc.numOFDMSymbolsPerFrame-2)*sc.numDataCarriers*log2(sc.M)-numBitsToSend;

% Add padding bits
bitsToSendWithPadding = [bitsToSend; randi([0,1], numBitsPadding, 1)];

% Transform bits into information symbols from M-QAM constellation
decInfoSymbols = bitsToDecSymbols(bitsToSendWithPadding, sc.M);
infoSymbols = modulator(decInfoSymbols, qammap(sc.M));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGNAL field building:
% It is a regular-size OFDM symbol used to transmit metadata. Here, it is
% used to transmit the number of useful bits contained in the frame (i.e.
% bits that are not padding)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Signal field uses BPSK modulation(+/-1)
numBitsToSendBin = de2bi(numBitsToSend, sc.numBitsForPayloadSize).';
signalSymbols = [numBitsToSendBin;...
    randi([0,1], sc.numDataCarriers-sc.numBitsForPayloadSize, 1)]*(-2)+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data frame building
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reshape data into matrix form ("-2" because of pilot and signal OFDM
% symbols)
dataFrame = reshape(infoSymbols, sc.numDataCarriers, sc.numOFDMSymbolsPerFrame-2);

% Pilot OFDM symbols (excluding pilot subcarriers) for channel estimation
% Randomization is needed to avoid high PAPR
pilotOfdmSymbol = buildPilotOfdmSymbol(sc.numUsedCarriers, sc.map, sc.M);

% Pilot subcarriers for CFO tracking and SFO correction
% Note: the first OFDM symbol correspond to the pilot OFDM symbol
pilotSubcarrier1 = repmat(pilotOfdmSymbol(1), 1, sc.numOFDMSymbolsPerFrame-1);
pilotSubcarrier2 = repmat(pilotOfdmSymbol(end), 1, sc.numOFDMSymbolsPerFrame-1);

% Add SIGNAL OFDM symbol to dataframe (it is considered as a regular OFDM
% symbol)
dataFrame = [signalSymbols(:), dataFrame];

% Add pilot subcarriers (for post-FFT frequency offset correction due to
% sampling frequency offset)
dataFrame = [pilotSubcarrier1;...
            dataFrame; ...
            pilotSubcarrier2];

% Add pilot OFDM symbol (for channel estimation)
dataFrame = [pilotOfdmSymbol, dataFrame];

% Add guard bands (unsued subcarriers)
% Outer subcarriers are severly attenuated due to SFO. 
% Not using them also makes guard bands to avoid frequency leakage
dataFrame = [zeros(sc.numZerosTop, sc.numOFDMSymbolsPerFrame); ...
            dataFrame(1:sc.numTotalCarriers/2-sc.numZerosTop,:); ...
            zeros(sc.zeroFreqSubcarrier, sc.numOFDMSymbolsPerFrame); ...
            dataFrame(sc.numTotalCarriers/2-sc.numZerosTop+1:end,:); ...
            zeros(sc.numZerosBottom, sc.numOFDMSymbolsPerFrame)];

% IFFT: fftshift() is needed to account for the negative frequencies since 
% the MATLAB function ifft() takes k=0,...,N-1 and not k=-N/2,...,0,...,N/2-1
dataFrameIFFT = ifft(fftshift(dataFrame,1), sc.NFFT);

% Add cyclic prefix (CP) to every OFDM symbol
dataFrameIFFTWithCp = [dataFrameIFFT(end-sc.CPLength+1:end,:);...
                  dataFrameIFFT];
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parallel to serial and concatenation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signalTx = [ca.';
           tenSts;...
           twoLtsWithCp;...
           dataFrameIFFTWithCp(:)];

end