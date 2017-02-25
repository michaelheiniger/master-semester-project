function [ fineDopplerEstimate ] = dopplerFineEstimate(signal, referenceSignal, Ts, coarseDopplerEstimate, fcorr)
%DOPPLERFINEESTIMATE Summary of this function goes here
%   Detailed explanation goes here

signalTruncated = signal(1:length(referenceSignal));

dopplerRange = [coarseDopplerEstimate-fcorr, coarseDopplerEstimate, coarseDopplerEstimate+fcorr]; 
innerProducts = zeros(length(dopplerRange), 1);

t = (0:(length(referenceSignal)-1))*Ts;
for k = 1:length(dopplerRange)
    signalCorrected = signalTruncated .* exp(-1j*2*pi*dopplerRange(k)*t);
    innerProducts(k) = sum(referenceSignal .* conj(signalCorrected));
end

d1 = coarseDopplerEstimate-fcorr;
d2 = coarseDopplerEstimate;
d3 = coarseDopplerEstimate+fcorr;

D = [ d1.^2, d1, 1 ; 
   d2.^2, d2, 1; 
   d3.^2, d3, 1 ];

f = [ innerProducts(1); innerProducts(2) ;innerProducts(3) ]; 
result = D\abs(f);
a = result(1);
b = result(2);

fineDopplerEstimate = -b/(2*a);