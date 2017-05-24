function decoderConfig()
%CONFIG Summary of this function goes here
%   Detailed explanation goes here

global decoder;

if ~isempty(decoder)
    return;
end

% Timing offset
% add error to timing estimate, in number of samples. 0 means that the
% estimate is untouched
decoder.timingOffset = 0;

% 1 if timing offset should be corrected
decoder.timingOffsetCorrection = 0;

% 1 if carrier frequency offset should be corrected
decoder.cfoCorrection = 1;

% 1 if residual CFO should be tracked and corrected
decoder.cfoTracking = 1;

% 1 if sampling frequency offset should be corrected
decoder.sfoCorrection = 0;

% 1 if channel equalization should be performed
decoder.equalization = 1;

% 1 if the received signal should be upsample for timing synchronization
decoder.upsample = 0;
decoder.USF = 10;

end

