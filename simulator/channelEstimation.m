function lambdas = channelEstimation(ofdmSymbolsRx, pilotOfdmSymbol, numUsedCarriers, numTotalCarriers, numZerosBottom, numZerosTop)

% Adapted from "channel_est3.m" file from SDR course taught by Prof. Rimoldi at EPFL
% CHANNEL_ESTIMATION Estimate the channel coefficients in the frequency domain
% There is no channel information available. Ky is estimated from the
% received data

% OFDMSYMBOLSRX: matrix containing the received data part of the frame
% PILOTOFDMSYMBOL1: the first training sequence
% NUM_CARRIERS: number of subcarriers (FFT/IFFT size)
% NUM_ZEROS: number of zero carriers at either end of the spectrum 

% LAMBDA: Column vector containing channel coefficients in the frequency domain

if numZerosBottom+numZerosTop == 0
   error('channelEstimation: There is no zero subcarriers: cannot measure the noise variance using MMSE !'); 
end

% Get the received pilot OFDM symbol
pilotOfdmSymbolRx = ofdmSymbolsRx(:,1);

% Measure noise variance using zero subcarriers
zeroSubcarriersTop = pilotOfdmSymbolRx(1:numZerosTop);
zeroSubcarriersBottom = pilotOfdmSymbolRx(end-numZerosBottom+1);
noiseVariance = 1/numTotalCarriers*var([zeroSubcarriersTop; zeroSubcarriersBottom]); % Normalization since noise is added in time domain

S = diag(pilotOfdmSymbol);
Kz = diag(noiseVariance*ones(numUsedCarriers,1)); % Diagonal since Z is assumed Gaussian independent

Ydata = ofdmSymbolsRx(:,2:end);

corMat = zeros(size(Ydata));
for i = 1:size(Ydata,2)
    Col = Ydata(:,i);
    xCol = xcorr(Col, 'biased');
    xCol = xCol(length(Col):end);
    corMat(:,i) = xCol;
end
Ky = toeplitz(mean(corMat,2));

lambdas = ((S\(Ky-Kz)) * (Ky\pilotOfdmSymbolRx));

% TO compare with the division
% lambda = ofdmSymbolsRx(:,1)./pilotOfdmSymbol;