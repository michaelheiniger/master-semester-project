function [coarseFrameRx, signalRx, frameBeginning] = rxTxUSRPFrameDetection(ofdmFramesTx, mode, Fs, usrpFrameLength, ofdmFrameLength, cpLength)
%RXTXUSRP Transmit and/or receive samples using USRP board(s). In receiving
% mode (see below), it can run indefinitely until it detects an OFDM 
% frame: the first OFDM frame detected in the received signal is extracted 
% and returned.
%
% There are 4 modes of operation:
% - oneBoardRx: The board is only used to receive samples (typically two
% boards and two computers are used: one board per computer, one in mode
% 'oneBoardRx' and the other in mode 'oneBoardTx')
% - oneBoardTx: The board is only used to transmit samples
% - loopback: One single board is used to transmit and receive samples: its
% transmit and receive antennas must be connected by a cable
% - twoBoardsRxTx: two boards are connected to a single computer. One of
% the board transmit, the other receive
% Repetition of code between the modes is intended: it is supposed to
% execute faster than factorized code.
%
% Parameters:
% - signalTx: the signal to be sent (can be empty for receiver modes)
% - mode: mode of operation of the communication system
% - Fs: Sampling frequency. Note: the USRP master clock rate is computer
% from this value and the interpolation factor. The interpolation factor
% may need to be adapted depending on the chosen Fs in order for clock rate
% to be in the valid range for the board.
% - usrpFrameLength: number of samples per USRP frame.
% - ofdmFrameLength: number of samples per OFDM frame: it is supposed to be
% smaller than samplesPerUsrpFrame
% - cpLength: length in samples of the cyclic prefix used for regular OFDM
% symbols (which is assumed to be larger than the length of the channel
% impulse response)
%
% Returns:
% - coarseFrameRx: the first OFDM frame detected in the received signal. It
% contains a cyclic prefix length extra samples to account for possible
% multipath (cyclic prefix length is assumed larger than channel impulse
% response)
% - signalRx: received samples (up to some fixed length since the receiver
% can run for a variable amount of time
% - frameBeginning: the position of the detected OFDM frame in the received
% signal (used for debugging typically)

if ofdmFrameLength > usrpFrameLength
    error('Length of OFDM frame must be smaller or equal than the length of USRP frame');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Board configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Serial number of the receiver board
serialNumberRx = '30C51BC';

% Serial number of the transmitter board
serialNumberTx = '30C5426';

% Sampling period
Ts = 1/Fs;

boardPlatform = 'B200'; % Model of the board ('B200' for B200 mini)
fc = 5e9; % carrier frequency
LOO = 1e3; % local oscillator offset
interpolationTx = 10; % Interpolation factor from host signal to USRP signal (e.g Fs = 0.5 Mhz => ClockRate = 5 Mhz)
clockRateTx  = interpolationTx/Ts; % main clock (Sampling rate of the digital signal sent to ADC) (5e6 to 56e6)
clockRateRx = clockRateTx; % on the same radio we need to use the same
decimationRx = interpolationTx; % to get to 1MHz
clockInputSource = 'Internal';
outputDataTypeUSRP = 'double'; % Difference between transport and output data type ?
gainTx = 89; %60 loopback 30dB attenuator; % 89dB over the air, no attenuator
gainRx = 50; % Max 76dB
% agc = comm.AGC('AdaptationStepSize', 0.01, 'AveragingLength', 10);
agc = comm.AGC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check connection with board(s)
% The goal here is to check only one every needed board depending on the
% mode since this process takes some time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
devicesToCheck = [];
if strcmp(mode, 'loopback')
    serialNumberRx = serialNumberTx;
    if usrpFrameLength > 1e5
        warning('samplePerFrame > 1e5 leads to overrun: the board drops samples.')
    end
    devicesToCheck = serialNumberRx;
elseif strcmp(mode, 'oneBoardTx')
    devicesToCheck = serialNumberTx;
elseif strcmp(mode, 'oneBoardRx')
    devicesToCheck = serialNumberRx;
elseif strcmp(mode, 'twoBoardsRxTx')
    devicesToCheck = [serialNumberRx; serialNumberTx];
end
% Actually check the connection with the device(s)
% checkDevicesConnection(devicesToCheck);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit and/or Receive using the board(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pad signal to be transmitted with zeros to fill the USRP frames
ofdmFramesTx = padWithZeros(ofdmFramesTx, usrpFrameLength);

% Compute the number of USRP frames needed to transmit the signal
noFrames = length(ofdmFramesTx)/usrpFrameLength
if mod(noFrames, 1) ~= 0
    error('Number of frames must be integer %s.', num2str(noFrames));
end

if strcmp(mode, 'loopback') || strcmp(mode, 'oneBoardTx') ||  strcmp(mode, 'twoBoardsRxTx')
    
    % Instantiate the Transmitter
    tx = comm.SDRuTransmitter(...
        'Platform', boardPlatform, ...
        'CenterFrequency', fc, ...
        'LocalOscillatorOffset', LOO,...
        'InterpolationFactor', interpolationTx, ...
        'MasterClockRate', clockRateTx,...
        'SerialNum', serialNumberTx,...
        'Gain', gainTx,...
        'ClockSource', clockInputSource,...
        'UnderrunOutputPort', true);
end

% Will contain the first N samples of the received signal (since the receiver
% waits for an OFDM frame, it can last indefinitely and thus a limit on
% the number of received samples saved is needed)
signalRx = zeros(2e6, 1);

% Will contain the first OFDM frame detected in the received signal
% The frame synchronization is coarse at this stage: due to multipath,
% the timing estimate should typically be at most a cyclic prefix
% length late
coarseFrameRx = zeros(ofdmFrameLength+cpLength, 1);

if strcmp(mode, 'loopback') || strcmp(mode, 'oneBoardRx') ||  strcmp(mode, 'twoBoardsRxTx')
    
    % Instantiate the Receiver
    rx = comm.SDRuReceiver(...
        'Platform', boardPlatform, ...
        'CenterFrequency', fc, ...
        'LocalOscillatorOffset', LOO,...
        'MasterClockRate', clockRateRx,...
        'SerialNum', serialNumberRx,...
        'DecimationFactor', decimationRx,...
        'Gain', gainRx,...
        'SamplesPerFrame', usrpFrameLength,...
        'ClockSource', clockInputSource,...
        'OutputDataType', outputDataTypeUSRP,...
        'OverrunOutputPort', true);
end

% Buffer used to account for the fact that an OFDM frame can span 2 USRP
% frames
lastSamplesBuffer = zeros(ofdmFrameLength, 1);

% Takes value 1 when an OFDM frame is found
frameFound = 0;

% Takes value 1 when the OFDM frame found is complete: used to account for
% the fact that an OFDM frame can span 2 USRP frames
frameComplete = 0;

% Save the number of missing samples of an OFDM frame that spans 2 USRP
% frames (i.e. the number of samples of the OFDM frame that must be in the
% next USRP frame)
numMissingSamples = 0;

% Save the position of the OFDM frame in the received signal. This is used
% only for debugging
frameBeginning = 0;

% Save the number of samples received: this is used to compute
% frameBeginning
numSamplesReceived = 0;

disp('Transmit and/or Receive ...');
switch mode
    
    % This mode transmits all the samples it has to transmit and then stops
    % whether an OFDM frame has been found or not. Since the received
    % samples are the transmitted ones, it does not make sense to wait
    % longer when all the samples have been transmitted
    case {'loopback', 'twoBoardsRxTx'}
        l = 1;
        while l <= noFrames
             
            % Receive samples for the receiver USRP
            [currentUsrpFrame, len, ~] = rx();
            
            if len > 0
                
                % Samples go through automatic gain controller: it is
                % needed since the frame detection algorithm uses an
                % absolute threshold
                numSamplesReceived = numSamplesReceived + length(currentUsrpFrame);
                
                currentUsrpFrame = agc(currentUsrpFrame);
                
                % Save received samples (for debugging)
                if l*usrpFrameLength < length(signalRx) 
                    signalRx(1+(l-1)*usrpFrameLength:l*usrpFrameLength) = currentUsrpFrame;
                end
                
                if not(frameFound) || not(frameComplete)
                    
                    % Try to detect an OFDM frame in the received samples
                    % (the last samples of the previous USRP frame are used
                    % in case the OFDM frame spans the previous and current
                    % USRP frames).
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
                        %disp(['num missing samples (frame was found): ' num2str(numMissingSamples)]); % for debugging
                        
                        % Index of first sample of the frame for plotting
                        % purpose only
                        frameBeginning = numSamplesReceived - length(currentUsrpFrame) - length(lastSamplesBuffer) + positionFromBeginning
                    end
                elseif frameFound && not(frameComplete)
                    
                    % Adds the second part of the detected OFDM frame
                    coarseFrameRx = [coarseFrameRx; currentUsrpFrame(1:numMissingSamples)];
                    
                    % The frame is now complete (it can span only 2 USRP
                    % frames)
                    frameComplete = 1; 
                end
                
                 % Transmit the samples to be transmitted
                tx(ofdmFramesTx(1+(l-1)*usrpFrameLength:l*usrpFrameLength));
                
                l = l+1;
                
                % Update buffer
                lastSamplesBuffer = currentUsrpFrame(end-length(lastSamplesBuffer)+1:end);
            end
        end
       
    % This mode looks for an OFDM frame in the received samples and keep
    % doing so until it find one (potentially indefinitely).
    case 'oneBoardRx'
        l = 1;
        while not(frameComplete)
             
            % Receive samples for the receiver USRP
            [currentUsrpFrame, len, ~] = rx();
            
            if len > 0
                
                % Samples go through automatic gain controller: it is
                % needed since the frame detection algorithm uses an
                % absolute threshold
                numSamplesReceived = numSamplesReceived + length(currentUsrpFrame);
                
                currentUsrpFrame = agc(currentUsrpFrame);
                
                % Save received samples (for debugging)
                if l*usrpFrameLength < length(signalRx) 
                    signalRx(1+(l-1)*usrpFrameLength:l*usrpFrameLength) = currentUsrpFrame;
                end
                
                if not(frameFound)
                    
                    % Try to detect an OFDM frame in the received samples
                    % (the last samples of the previous USRP frame are used
                    % in case the OFDM frame spans the previous and current
                    % USRP frames).
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
                        %disp(['num missing samples (frame was found): ' num2str(numMissingSamples)]); % for debugging
                        
                        % Index of first sample of the frame for plotting
                        % purpose only
                        frameBeginning = numSamplesReceived - length(currentUsrpFrame) - length(lastSamplesBuffer) + positionFromBeginning
                    end
                elseif frameFound && not(frameComplete)
                    
                    % Adds the second part of the detected OFDM frame
                    coarseFrameRx = [coarseFrameRx; currentUsrpFrame(1:numMissingSamples)];
                    
                    % The frame is now complete (it can span only 2 USRP
                    % frames)
                    frameComplete = 1; 
                end
                
                l = l+1;
                
                % Update buffer
                lastSamplesBuffer = currentUsrpFrame(end-length(lastSamplesBuffer)+1:end);
            end
        end
        
    % Transmits the samples to transmit and stops
    case 'oneBoardTx'
        for l = 1:noFrames
            tx(ofdmFramesTx(1+(l-1)*usrpFrameLength:l*usrpFrameLength));
            coarseFrameRx = [];
        end
    otherwise
        error('Unknown mode of operation.')
end

disp('End of transmission and / or reception');

% Release transmitter ressources
if strcmp(mode, 'loopback') || strcmp(mode, 'oneBoardTx') ||  strcmp(mode, 'twoBoardsRxTx')
    release(tx);
end

% Release receiver ressources
if strcmp(mode, 'loopback') || strcmp(mode, 'oneBoardRx') ||  strcmp(mode, 'twoBoardsRxTx')
    release(rx);
end   

end

