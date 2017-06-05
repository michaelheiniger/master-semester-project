function h = createMultipathChannelFilter(amplitudes, delays)
% Adapted from "create_multipath_channel_filter.m" file from SDR course taught by Prof. Rimoldi at EPFL

% CREATEMULTIPATHCHANNELFILTER(AMPLITUDES, DELAYS) Create sampled response of multipath channel
% H = CREATEMULTIPATHCHANNELFILTER(AMPLITUDES, DELAYS)
%
% We assume that the shaping pulse is square. In this case, the
% autocorrelation function will be a triangle.
%
% DELAYS and AMPLITUDES are vectors of the same length, specifying the
% strength and delay of each path. The DELAYS must be specified relative to
% the sampling period.
%
% H contains the samples of filtered impulse response.

if (numel(amplitudes) ~= numel(delays))
    error('createMultipathChannelFilter:wrongInputDimensions', 'AMPLITUDES and DELAYS must be vectors of the same length');
end
        
% this is the length of the channel in symbols (tmax integer or not)    
L = max(ceil(delays))+1; 

M = numel(delays);

% This implementation is the same as the alternative one proposed
% later, but it is just more compact. It uses the fact that a certain delay
% can affect only two consecutive Rp (Rp is a triangle in our case).
h = zeros(1,L);
for i = 1:M
    n = floor(delays(i)) + 1; % +1 because we cannot index h at 0 (Matlab stuff :-) )
    % Draw the two triangles at n and (n+1) to see where is n and how we compute the
    % values below
    h(n) = h(n) + amplitudes(i)*(n-delays(i));
    if n < length(h) % takes care of the case where tauMax is integer (otherwise we introduce a non-necessary 0 at the end)
        h(n+1) = h(n+1) + amplitudes(i)*(delays(i)-n+1);
    end
end

end