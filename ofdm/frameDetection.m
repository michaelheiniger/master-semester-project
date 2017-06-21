function [coarseFrame, frameFound, frameComplete, positionFromBeginning] = frameDetection(currentChunk, ofdmFrameLength, cpLength)
%FRAMEDETECTION Tries to detect the next OFDM frame encountered in the
%received signal: stops when found
% Params:
% - currentChunk: current USRP frame with last samples from the previous
% USRP frame
% - ofdmFrameLength: the length of a complete OFDM frame
% - cpLength: the length of the cyclic prefix (which is assumed larger than
% the channel impulse response)
%
% Returns:
% - coarseFrame: the OFDM frame detected with the "cyclic prefix length"
% last samples to account for multipath (i.e. this is a coarse frame sync)
% - frameFound: returns 1 if an OFDM frame was found and 0 otherwise.
% - frameComplete: returns 1 if a complete OFDM frame was found (i.e. no
% samples of the OFDM frame is contained in the next USRP frame)
% - positionBeginning: position of the detected OFDM frame in currentChunk
% (for debugging)

[stsTime, ~] = getSTS();
tenStsTime = repmat(stsTime, 10, 1);

L = length(tenStsTime);
K = ofdmFrameLength-L;

% Threshold used to decide if an OFDM frame is detected or not
threshold = 12;

% Save the detected OFDM frame
coarseFrame = zeros(ofdmFrameLength+cpLength, 1);

% 1 if OFDM frame is found
frameFound = 0;

% 1 if OFDM frame detected is complete (i.e. no samples are in the next
% USPR frame)
frameComplete = 0;

% Position of the OFDM frame from the beginning of the currentChunk
positionFromBeginning = 0;

% Position of the maximum value of the cross-correlation
peakPosition = 1;

% Loop needed to find the FIRST encountered OFDM frame
for l = 1:K:length(currentChunk)

    if length(currentChunk(l:end)) < K+L-1
        break;
    end
    
    currentSamples = currentChunk(l:l+K-1);
    
    [R1, ~] = xcorr(currentSamples, tenStsTime);
    R1 = R1(length(currentSamples):end);
    [val1, pos1] = max(abs(R1));

    if val1 > threshold
        disp('Threshold reached')
        
        % Frame is found, check if it is completely contained in
        % currentChunk, else the frame will be found next iteration
        % In both case, we need to break the loop.
        % +(L-1) is needed in case the 10STS (of length L) are not complete
        % we check if there is a higher peak when including the next (L-1)
        % samples. If there is, then it means that the 10STS were not 
        % completely in currentSamples

        % Take [currentSamples; (L-1) next samples]
        currentSamplesExtended = currentChunk(l:l+K-1+(L-1));
        [R2, ~] = xcorr(currentSamplesExtended, tenStsTime);
        R2 = R2(length(currentSamplesExtended):end);
        [val2, pos2] = max(abs(R2));

        frameFound = 1;
        disp('Frame detected')

        if val1 > val2
            peakPosition = l+pos1-1;
        else
            peakPosition = l+pos2-1;
        end

        plotSignalMagnitude(R1,'samples','xcorr with first part');
        plotSignalMagnitude(R2,'samples','xcorr with second part');

        % If the OFDM frame is contained entirely in this current USRP
        % frame
        if length(currentChunk(l+peakPosition-1:end)) >= ofdmFrameLength
            % -cpLength to account for CIR length of at most cpLength
            coarseFrame = currentChunk(peakPosition-cpLength:peakPosition+ofdmFrameLength-1);
            length(coarseFrame)
            frameComplete = 1;
        else % Else, we take all remaining samples and we will need to read the next USRP frame
            coarseFrame = currentChunk(l+peakPosition-cpLength:end);
            frameComplete = 0;
        end
        positionFromBeginning = peakPosition;

        % We want only the first encountered frame
        break;
    end
end

end