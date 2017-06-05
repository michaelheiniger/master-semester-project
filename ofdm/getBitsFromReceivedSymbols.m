function [bitsReceived] = getBitsFromReceivedSymbols(infoSymbolsRx, map, M)
%GETBITSFROMRECEIVEDSYMBOLS Demodulates the symbols into bits

% Demodulates symbols estimate (Hard Decision)
decInfoSymbolsEst = demodulator(infoSymbolsRx, map);

% Convert into bits
bitsReceived = decSymbolsToBits(decInfoSymbolsEst, M);

end

