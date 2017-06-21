function [coarseFrameRx, signalRx, frameBeginning] = usrpSimulator(ofdmFramesTx, usrpFrameLength, ofdmFrameLength, cpLength, deterministic)
%USRPSIMULATOR Simulates the behavior of transmission/reception using USRPs.
% OFDM frames are randomly embedded in USRP frames: an OFDM frame may span
% 2 USRP frames, not more. Indeed, the number of samples per USRP frame 
% (i.e. samplesPerFrames) is assumed to be bigger than the length of one 
% OFDM frame.
% Simulates the reception of the signal and extract the first OFDM frame
% encountered in the signal (coarse frame sync). The frame is returned 
% with the previous CpLength samples (to account for ISI) for 
% fine frame sync
% 
% Params:
% - ofdmFramesTx: Concatenation of OFDM frames to send
% - usrpFrameLength: number of samples per USRP frame
% - ofdmFrameLength: number of samples per OFDM frame
% - cpLength: length in samples of the cyclic prefix of regular OFDM
% symbols
% - deterministic: 1: the OFDM frames are placed at the beginning of the
% received USRP frames. 0: the OFDM frames are placed at a random position
% in the USRP frames, that is, the first sample of the first OFDM frame is
% random, all other OFDM frames follows directly.

if ofdmFrameLength > usrpFrameLength
    error('Length of OFDM frame must be smaller or equal than the length of USRP frame');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build and fill USRP frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numUsrpFramesNeeded = ceil(usrpFrameLength/length(ofdmFramesTx));

% Add randomness in the number of USRP frames received
numUsrpFramesUsed = numUsrpFramesNeeded + randi([0 10]);

totalNumSamplesRx = numUsrpFramesUsed * usrpFrameLength;

signalRx = zeros(totalNumSamplesRx, 1);

if deterministic
    signalRx(1:length(ofdmFramesTx)) = ofdmFramesTx;
else
    % Add randomness in where the OFDM frames are located within USRP frames
    frameBeginningMaxValue = totalNumSamplesRx - length(ofdmFramesTx) + 1
    frameBeginning = randi([1 frameBeginningMaxValue], 1, 1);
    signalRx(frameBeginning:frameBeginning+length(ofdmFramesTx)-1) = ofdmFramesTx;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulates the reception of the signal and extract the first OFDM frame
% encountered in the signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coarseFrameRx = zeros(ofdmFrameLength+cpLength, 1);
lastSamplesBuffer = zeros(ofdmFrameLength, 1);
frameFound = 0;
frameComplete = 0;
numMissingSamples = 0;
frameBeginning = 0;    
numSamplesReceived = 0;
    
lengthSignalRx = length(signalRx);

agc = comm.AGC;

l = 1;
while not(frameComplete)

    if l*usrpFrameLength > lengthSignalRx
        error('Received signal does not contain a complete OFDM frame');
    end
            
    currentUsrpFrame = signalRx(1+(l-1)*usrpFrameLength : l*usrpFrameLength);
    numSamplesReceived = numSamplesReceived + length(currentUsrpFrame);
    
    currentUsrpFrame = agc(currentUsrpFrame);

    if not(frameFound)
        % frameComplete is 1 if the OFDM frame is contained in
        % a single USRP frame, else, the next USRP frame needs
        % to be read to extract the missing samples
        [coarseFrame, frameFound, frameComplete, positionFromBeginning] = frameDetection([lastSamplesBuffer; currentUsrpFrame], ofdmFrameLength, cpLength);

        if frameFound
            % Save the first part of the OFDM frame received
            % (possibly the whole OFDM frame)
            coarseFrameRx(1:length(coarseFrame)) = coarseFrame;

            % Save the number of missing samples
            numMissingSamples = length(coarseFrameRx) - length(coarseFrame);
            disp(['num missing samples (frame was found): ' num2str(numMissingSamples)]);

            % Index of first sample of the frame (for debugging)
            frameBeginning = numSamplesReceived - length(currentUsrpFrame) - length(lastSamplesBuffer) + positionFromBeginning;
        end
    elseif frameFound && not(frameComplete)

        coarseFrameRx = [coarseFrameRx; currentUsrpFrame(1:numMissingSamples)];
        frameComplete = 1;
    end
    l = l+1;

    % Update buffer
    lastSamplesBuffer = currentUsrpFrame(end-length(lastSamplesBuffer)+1:end);
end

end

