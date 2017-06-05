function [bits] = textFileToBits(filePath)
%TEXTFILETOBITS read text file and returns the corresponding bits (wrt
% ASCII)

fid = fopen(filePath);

% Read characters
chars = fread(fid,'*char');

% Convert to bits
charInBits = transpose(dec2bin(chars,8));

% Serialization
bits = charInBits(:)-'0'; % -'0' is needed to change format from char to int

end

