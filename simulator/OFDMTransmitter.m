function [signalTx, dataFrame] = OFDMTransmitter(config, mode, dataSymbolsTx, stsTime, ltsTime, pilotOfdmSymbol, pilots1, pilots2)

c = config;

% Build 10 STS
tenSts = repmat(stsTime, 10, 1);

% Build 2 LTS
twoLts = repmat(ltsTime, 2, 1);

% Add cyclic prefix (CP) to the two LTS
twoLtsWithCp = [twoLts(end-c.twoLtsCpLength+1:end); twoLts];

% OFDM data frame construction
dataFrame = reshape(dataSymbolsTx, c.numDataCarriers, c.numOFDMSymbolsPerFrame-1);

% Add pilot OFDM symbol (for channel estimation)
dataFrame = [pilotOfdmSymbol, dataFrame];

% Add pilot subcarriers (for post-FFT frequency offset correction due to
% sampling frequency offset)
dataFrame = [pilots1;...
            dataFrame; ...
            pilots2];

% Add guard bands (unsued subcarriers)
% Outer subcarriers are severly attenuated due to SFO. 
% Not using them also makes guard bands to avoid frequency leakage
dataFrame = [zeros(c.numZerosTop, c.numOFDMSymbolsPerFrame); ...
            dataFrame(1:c.numTotalCarriers/2-c.numZerosTop,:); ...
            zeros(c.dcSubcarrier, c.numOFDMSymbolsPerFrame); ...
            dataFrame(c.numTotalCarriers/2-c.numZerosTop+1:end,:); ...
            zeros(c.numZerosBottom, c.numOFDMSymbolsPerFrame)];

% IFFT: fftshift() is needed to account for the negative frequencies since 
% the MATLAB function ifft() takes k=0,...,N-1 and not k=-N/2,...,0,...,N/2-1
dataFrameIFFT = ifft(fftshift(dataFrame,1), c.NFFT);

% Add cyclic prefix (CP) to every OFDM symbol
dataFrameIFFTWithCp = [dataFrameIFFT(end-c.CPLength+1:end,:);...
                  dataFrameIFFT];

% Serialization
signalTx = [tenSts;...
           twoLtsWithCp;...
           dataFrameIFFTWithCp(:)];
       
% Add garbage before and after OFDM frame and repeat
if strcmp(mode, 'simulation')
    signalTx = repmat([zeros(1000,1); signalTx; zeros(1000,1)],1,1);
else
%     signalTx = repmat([complex(0.001+1i*0.001+zeros(1.5e5,1)); signalTx; complex(0.001+1i*0.001+zeros(1.5e5,1))],10,1);
    signalTx = repmat(signalTx,1000,1);
end

end

