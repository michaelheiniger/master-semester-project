function lambda = channelEstimation(ofdmSymbolsRx, pilotOfdmSymbol, numTotalCarriers, numZeros)

% Adapted from "channel_est3.m" file from SDR course taught by Prof. Rimoldi at EPFL
% CHANNEL_ESTIMATION Estimate the channel coefficients in the frequency domain
% There is no channel information available. Ky is estimated from the
% received data

% OFDMSYMBOLSRX: matrix containing the received data part of the frame
% PILOTOFDMSYMBOL1: the first training sequence
% NUM_CARRIERS: number of subcarriers (FFT/IFFT size)
% NUM_ZEROS: number of zero carriers at either end of the spectrum 

% LAMBDA: Column vector containing channel coefficients in the frequency domain

pilotOfdmSymbolRx = ofdmSymbolsRx(:,1);
noiseVariance = var(pilotOfdmSymbol-pilotOfdmSymbolRx); %% IS THAT CORRECTLY ESTIMATED ???

S = diag(pilotOfdmSymbol);
Kz = diag(noiseVariance*ones(numTotalCarriers-2*numZeros,1));

Ydata = ofdmSymbolsRx(:,2:end-1);

corMat = zeros(size(Ydata));
for i = 1:size(Ydata,2)
    Col = Ydata(:,i);
    xCol = xcorr(Col,'biased');
    xCol = xCol(length(Col):end);
    corMat(:,i) = xCol;
end

Ky = toeplitz(mean(corMat,2));
lambda = ((S\(Ky-Kz)) * (Ky\pilotOfdmSymbolRx));
  
  


