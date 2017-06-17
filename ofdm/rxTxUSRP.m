function [signalRx] = rxTxUSRP(signalTx, mode, Fs)
%RXTXUSRP Summary of this function goes here
%   Detailed explanation goes here

[signalTx, wasRowVector] = turnIntoRowVector(signalTx);

serialNumberRx = '30C51BC';
serialNumberTx = '30C5426';

Ts = 1/Fs;

boardPlatform = 'B200';
fc = 2.49e9; % carrier frequency
% fc = 100e6; % carrier frequency
LOO = 1e3; % local oscillator offset
% LOO = 1e3; % local oscillator offset
interpolationTx = 10; % Interpolation factor from host signal to USRP signal (e.g Fs = 0.5 Mhz => ClockRate = 5 Mhz)
clockRateTx = interpolationTx/Ts; % main clock (Sampling rate of the digital signal sent to ADC) (5e6 to 56e6)
clockRateRx = clockRateTx; % on the same radio we need to use the same
decimationRx = interpolationTx; % to get to 1MHz
clockInputSource = 'Internal';
outputDataTypeUSRP = 'double'; % Difference between transport and output data type ?
gainTx = 89; %60 loopback 30dB attenuator; % 89 over the air, no attenuator
gainRx = 70; % Max 76dB

% Number of samples per USRP frame
samplesPerFrame = 0.5e5;

burstMode = false;

% Ts = interpolationTx/clockRateTx; % [s] sampling time
Fs = 1/Ts;

% Pad signal with zeros
signalTx = padWithZeros(signalTx, samplesPerFrame);

% Compute the number of USRP frames needed to transmit the signal
noFrames = length(signalTx)/samplesPerFrame;
if mod(noFrames, 1) ~= 0
    error('Number of frames must be integer %s.', num2str(noFrames));
end

devicesToCheck = [];
if strcmp(mode, 'loopback')
    serialNumberRx = serialNumberTx;
    if samplesPerFrame > 1e5
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

% checkDevicesConnection(devicesToCheck);

% Instantiate transmitter object
if strcmp(mode, 'loopback') || strcmp(mode, 'oneBoardTx') ||  strcmp(mode, 'twoBoardsRxTx')
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

% Instantiate receiver object
if strcmp(mode, 'loopback') || strcmp(mode, 'oneBoardRx') ||  strcmp(mode, 'twoBoardsRxTx')
    
    signalRx = zeros(1, noFrames*samplesPerFrame);
    
    rx = comm.SDRuReceiver(...
        'Platform', boardPlatform, ...
        'CenterFrequency', fc, ...
        'LocalOscillatorOffset', LOO,...
        'MasterClockRate', clockRateRx,...
        'SerialNum', serialNumberRx,...
        'DecimationFactor', decimationRx,...
        'Gain', gainRx,...
        'SamplesPerFrame', samplesPerFrame,...
        'ClockSource', clockInputSource,...
        'OutputDataType', outputDataTypeUSRP,...
        'OverrunOutputPort', true);
end

disp('Transmit and/or Receive ...');
switch mode
    case 'twoBoardsRxTx'
        l = 1;
        while l <= noFrames
            [dataRxUSRP, len, ~] = rx();
            
            if len > 0
                signalRx(1+(l-1)*samplesPerFrame:l*samplesPerFrame) = dataRxUSRP.';
                tx(signalTx(1+(l-1)*samplesPerFrame:l*samplesPerFrame).');
                l = l+1;
            end
        end
        
    case 'loopback' % also for 'twoBoardsRxTx'
        l = 1;
        while l <= noFrames
            [dataRxUSRP, len, ~] = rx();
            
            if len > 0
                signalRx(1+(l-1)*samplesPerFrame:l*samplesPerFrame) = dataRxUSRP.';
                tx(signalTx(1+(l-1)*samplesPerFrame:l*samplesPerFrame).');
                l = l+1;
            end
        end
        
    case 'oneBoardRx'
        l = 1;
        while l <= 160
            [dataRxUSRP, len, ~] = rx();
            
            if len > 0
                signalRx(1+(l-1)*samplesPerFrame:l*samplesPerFrame) = dataRxUSRP.';
                l = l+1;
            end
        end
    case 'oneBoardTx'
        for l = 1:noFrames
            tx(signalTx(1+(l-1)*samplesPerFrame:l*samplesPerFrame).');
            signalRx = [];
        end
    otherwise
        error('Unknown mode of operation.')
end

signalRx = changeVectorBackToColumnIfNeeded(signalRx, wasRowVector);

disp('End of transmission and / or reception');

if strcmp(mode, 'loopback') || strcmp(mode, 'oneBoardTx') ||  strcmp(mode, 'twoBoardsRxTx')
    release(tx);
end

if strcmp(mode, 'loopback') || strcmp(mode, 'oneBoardRx') ||  strcmp(mode, 'twoBoardsRxTx')
    release(rx);
end   

end

