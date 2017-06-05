% dataFrame is returned for experiment purposes only (plots, ...)
function [signalTx, dataFrame] = OFDMTransmitter(systemConfig, dataSymbolsTx, stsTime, ltsTime, signalSymbols, pilotOfdmSymbol, pilots1, pilots2, ca)

sc = systemConfig;

% Build 10 STSs vector
tenSts = repmat(stsTime, 10, 1);

% Build 2 LTSs vector
twoLts = repmat(ltsTime, 2, 1);

% Add cyclic prefix to the two LTS
twoLtsWithCp = [twoLts(end-sc.twoLtsCpLength+1:end); twoLts];

% Build SIGNAL OFDM symbol
signalOfdmSymbol = [zeros(sc.numZerosTop, 1); ...
                    signalSymbols(1:sc.numTotalCarriers/2-sc.numZerosTop); ...
                    zeros(sc.zeroFreqSubcarrier, 1); ...
                    signalSymbols(sc.numTotalCarriers/2-sc.numZerosTop+1:end); ...
                    zeros(sc.numZerosBottom, 1)];
% Apply IFFT
signalOfdmSymbolIfft = ifft(fftshift(signalOfdmSymbol), 64);

% Add CP
signalOfdmSymbolIfftWithCp = [signalOfdmSymbolIfft(end-sc.signalFieldCpLength+1:end); signalOfdmSymbolIfft];

% OFDM data frame construction
dataFrame = reshape(dataSymbolsTx, sc.numDataCarriers, sc.numOFDMSymbolsPerFrame-1);

% Add pilot subcarriers (for post-FFT frequency offset correction due to
% sampling frequency offset)
dataFrame = [pilots1;...
            dataFrame; ...
            pilots2];

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
              
% Serialization
signalTx = [ca.';
           tenSts;...
           twoLtsWithCp;...
           signalOfdmSymbolIfftWithCp;
           dataFrameIFFTWithCp(:)];

end