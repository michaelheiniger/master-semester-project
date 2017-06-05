function lambdas = channelEstimation(ofdmSymbolsRx, pilotOfdmSymbolRx, pilotOfdmSymbol, numUsedCarriers, noiseVariance)
% Adapted from "channel_est3.m" file from SDR course taught by Prof. Rimoldi at EPFL

% CHANNELESTIMATION Estimate the channel coefficients in the frequency domain
% There is no channel information available. Ky is estimated from the
% received data
% Inputs:
% - ofdmSymbolsRx: matrix containing the received data part of the frame
% - pilotOfdmSymbolRx: the received pilot OFDM symbol WITH guard bands
% - pilotOfdmSymbol: the clean pilot OFDM symbol
% - numUsedCarriers: number of subcarriers used to carry symbols (information
%   or pilot)
% - noiseVariance: noise variance estimate

S = diag(pilotOfdmSymbol);
Kz = diag(noiseVariance*ones(numUsedCarriers,1)); % Diagonal since Z is assumed Gaussian independent

Ydata = ofdmSymbolsRx;

corMat = zeros(size(Ydata));
for i = 1:size(Ydata,2)
    Col = Ydata(:,i);
    xCol = xcorr(Col, 'biased');
    xCol = xCol(length(Col):end);
    corMat(:,i) = xCol;
end
Ky = toeplitz(mean(corMat,2));

lambdas = ((S\(Ky-Kz)) * (Ky\pilotOfdmSymbolRx));