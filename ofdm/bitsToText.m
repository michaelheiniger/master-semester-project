function [text] = bitsToText(bits)
%BITSTOTEXT return the text corresponding to the provided bits

if mod(length(bits), 8) ~= 0
   numCharsToKeep = floor(length(bits)/8);
   bits = bits(1:numCharsToKeep*8);
end

bitsGrouped = transpose(reshape(bits, 8, []));

text = transpose(char(bin2dec(num2str(bitsGrouped))));

end

