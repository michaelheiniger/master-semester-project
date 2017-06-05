function [bits] = decSymbolsToBits(decSymbols, M)
%DECSYMBOLSTOBITS Summary of this function goes here
%   Detailed explanation goes here

decSymbols = decSymbols(:);

bits = transpose(de2bi(decSymbols));

bits = bits(:);

end

