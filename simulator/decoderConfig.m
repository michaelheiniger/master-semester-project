function decoderConfig()
%CONFIG Summary of this function goes here
%   Detailed explanation goes here

global decoder;

if ~isempty(decoder)
    return;
end

% 1 if frequency correction should be done
decoder.frequencyCorrection = 1;

% 1 if the received signal should be upsample for timing synchronization
decoder.upsample = 0;
decoder.USF = 10;

end

