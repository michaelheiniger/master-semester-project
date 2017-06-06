function [ofdmSymbolsCorrected] = carrierFrequencyOffsetTracking(ofdmSymbolsRx, pilots1, pilots2)
%CARRIERFREQUENCYOFFSETTRACKING This function uses the pilot symbols
%present in each OFDM symbol to track and remove the residual CFO

pilots1Rx = ofdmSymbolsRx(1,:);
pilots2Rx = ofdmSymbolsRx(end,:);

% Compute the phase of pilot symbols for every OFDM symbols
phase1 = angle(pilots1Rx.*conj(pilots1));
phase2 = angle(pilots2Rx.*conj(pilots2));

% Compute the mean of the phase of the two pilots for every OFDM symbols
meanPhase = (phase1+phase2)/2;

% Remove the phase from all OFDM symbols
ofdmSymbolsCorrected = ofdmSymbolsRx .* exp(-1j*meanPhase);

end

