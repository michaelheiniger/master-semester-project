function [ofdmSymbolsCorrected] = carrierFrequencyOffsetTracking(ofdmSymbols, pilots1, pilots2)
%CARRIERFREQUENCYOFFSETTRACKING This function uses the pilot symbols
%present in each OFDM symbol to track and remove the residual CFO

pilots1Rx = ofdmSymbols(1,:);
pilots2Rx = ofdmSymbols(end,:);

pilotOfdmSymbolRx = ofdmSymbols(:,1);
pilots1 = [pilotOfdmSymbolRx(1), pilots1];
pilots2 = [pilotOfdmSymbolRx(end), pilots2];

% Compute the phase of pilot symbols for every OFDM symbols
phase1 = angle(pilots1Rx.*conj(pilots1));
phase2 = angle(pilots2Rx.*conj(pilots2));

% Compute the mean of the phase of the two pilots for every OFDM symbols
meanPhase = (phase1+phase2)/2;

% Remove the phase from all OFDM symbols
ofdmSymbolsCorrected = ofdmSymbols .* exp(-1j*meanPhase);

end

