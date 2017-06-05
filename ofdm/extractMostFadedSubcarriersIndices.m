function [lowestChannelsCoeffsIndices] = extractMostFadedSubcarriersIndices(channelCoefficients, threshold, numUsedCarriers)
%EXTRACTMOSTFADEDSUBCARRIERSINDICES Summary of this function goes here
%   Detailed explanation goes here

lowestChannelCoeffs = quantile(abs(channelCoefficients), [threshold, 1]);
disp([num2str(threshold) '% of subcarriers have a magnitude below ' num2str(lowestChannelCoeffs(1))]);
lowestChannelsCoeffsIndices = find(abs(channelCoefficients) <= lowestChannelCoeffs(1)).';

negativeIndices = lowestChannelsCoeffsIndices(lowestChannelsCoeffsIndices < numUsedCarriers/2)-numUsedCarriers/2-1;
positiveIndices = lowestChannelsCoeffsIndices(lowestChannelsCoeffsIndices >= numUsedCarriers/2)-numUsedCarriers/2;
disp(['Most faded subcarriers (' num2str(threshold) '): ' num2str([negativeIndices, positiveIndices])]);

end

