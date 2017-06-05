function [BER, SER] = symbolsAndBitsStats(infoSymbols, infoSymbolsRx, bitsSent, bitsReceived, map)
%SYMBOLSANDBITSSTATS Summary of this function goes here
%   Detailed explanation goes here

% Demodulates symbols estimate (Hard Decision)
decInfoSymbolsEst = demodulator(infoSymbolsRx, map);

% Convert received info symbols into their M-QAM symbols estimate
infoSymbolsEst = modulator(decInfoSymbolsEst, map);

% Compute Symbol Error Rate (SER)
SER = sum(infoSymbols ~= infoSymbolsEst)/length(infoSymbols);

% Compute Bit Error Rate (BER)
BER = sum(bitsSent ~= bitsReceived)/length(bitsSent);

end

