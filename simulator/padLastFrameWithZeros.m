function [frame] = padLastFrameWithZeros(frame, samplesPerFrame)
%PADLASTFRAMEWITHZEROS Summary of this function goes here
%   Detailed explanation goes here

paddingLength = ceil(length(frame)/samplesPerFrame)*samplesPerFrame-length(frame);
frame = [frame zeros(1, paddingLength)];

end

