function [decSymbols] = bitsToDecSymbols(bits, M)
%BITSTODECSYMBOLS turns vectors of bits into symbols as decimals {0,...,M-1}

bits = bits(:);

bitsPerSymbol = log2(M);

bitsGrouped = transpose(reshape(bits, bitsPerSymbol, []));

decSymbols = bi2de(bitsGrouped);

end

