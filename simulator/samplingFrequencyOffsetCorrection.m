function [ofdmSymbolsCorrected] = samplingFrequencyOffsetCorrection(ofdmSymbols, pilots1, pilots2, numZerosTop, numZerosBottom, numTotalCarriers, dcSubcarrier)
%SAMPLINGFREQUENCYOFFSETESTIMATION Summary of this function goes here
%   Detailed explanation goes here

pilots1Rx = ofdmSymbols(1,:);
pilots2Rx = ofdmSymbols(end,:);

% Compute the phase of pilot symbols for every OFDM symbols
phase1 = angle(pilots1Rx.*conj(pilots1));
phase2 = angle(pilots2Rx.*conj(pilots2));

% Compute the slope of the phase line or every OFDM symbol
slopes = (phase2-phase1) / (numTotalCarriers-numZerosTop-numZerosBottom-1);

% Compute the x-axis indices for every subcarriers
if dcSubcarrier % DC subcarrier DOES NOT contain data
    offsetFromPhaseRef = [0:numTotalCarriers/2-numZerosTop-1, numTotalCarriers/2-numZerosTop+1:numTotalCarriers-numZerosTop-numZerosBottom-1];
else % DC subcarrier DOES contain data
    offsetFromPhaseRef = 0:(numTotalCarriers-numZerosTop-numZerosBottom-1);
end

% Repeat the x-axis indices for every OFDM symbols
K = repmat(offsetFromPhaseRef.', 1, size(pilots1Rx,2));

% Compute the phases for all subcarriers and all OFDM symbols
phases = phase1 + slopes.*K;

% Remove the phase from all OFDM symbols
ofdmSymbolsCorrected = ofdmSymbols .* exp(-1j*phases);

end

