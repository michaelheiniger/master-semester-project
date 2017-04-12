function [signalTx, dataFrame] = OFDMTransmitter(config, mode, symbolsTx, shortTrainingSeqIFFT, longTrainingSeqIFFT, pilotOfdmSymbol, pilots1, pilots2)

c = config;

% First preamble construction
preamble1 = repmat(shortTrainingSeqIFFT, 10, 1);

% Second preamble construction
preamble2 = repmat(longTrainingSeqIFFT, 2, 1);

% Add CP to second preamble
secondPreambleIFFTCP = [preamble2(end-c.preamble2CPLength+1:end); preamble2];

% OFDM data frame construction
dataFrame = reshape(symbolsTx, c.numDataCarriers-c.numPilots, c.numOFDMSymbolsPerFrame-1);

% Add pilot subcarriers (for post-FFT frequency offset correction due to
% sampling frequency offset)

dataFrame = [pilots1;...
            dataFrame; ...
            pilots2];

% Add pilot ofdm symbol (for channel estimation)
dataFrame = [pilotOfdmSymbol, dataFrame];

dataFrame = [zeros(c.numZeros, c.numOFDMSymbolsPerFrame); ...
            dataFrame; ...
            zeros(c.numZeros, c.numOFDMSymbolsPerFrame)];

% Apply ffshift to account for the negative frequencies since the MATLAB
% function ifft takes k=0,...,N-1 and not k=-N/2,...,0,...,N/2-1
% Apply ifft to get the time domain signal
dataFrameIFFT = ifft(fftshift(dataFrame,1), c.NFFT);
% Add CP to data frame
dataFrameIFFTCP = [dataFrameIFFT(end-c.CPLength+1:end,:);...
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



end

