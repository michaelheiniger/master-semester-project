function [data] = padWithZeros(data, totalLength)
%PADLASTFRAMEWITHZEROS Summary of this function goes here
%   Detailed explanation goes here

paddingLength = ceil(length(data)/totalLength)*totalLength-length(data);

[~, cols] = size(data);
if cols ~= 1
    padding = zeros(1, paddingLength);
    data = [data padding];
else
    padding = zeros(paddingLength, 1);
    data = [data; padding];
end

end

