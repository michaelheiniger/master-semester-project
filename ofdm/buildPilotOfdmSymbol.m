function [pilotOfdmSymbol] = buildPilotOfdmSymbol(numSymbols, map, M)
%BUILDPILOTOFDMSYMBOL Summary of this function goes here
%   Detailed explanation goes here

% The symbols used as pilot must be known by the transmitter and the
% receiver. Since the purpose of this setup is to perform experiments, one
% can change parameters that will influence the length of the pilot OFDM
% symbol. To avoid precomputing tables of pilot OFDM symbols for all the
% possible parameters, we can use a structure way of generating the pilot
% symbols. Here, we take all the symbols of the constellation and repeat
% them as many times as needed to fill the OFDM symbol. The receiver and
% transmitter don't need any prior communication (except having the same
% setup parameters) to generate the same pilot OFDM symbol.

symbols = map(repmat(1:M, 1, ceil(numSymbols/M)));

% Take as many symbols as needed
pilotOfdmSymbol = symbols(1:numSymbols);

pilotOfdmSymbol = pilotOfdmSymbol(:);
end

