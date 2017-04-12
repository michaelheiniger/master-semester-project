function [dataRx,Fs] = rxTxUSRP(dataTx, samplesPerFrame, mode)
%RXTXUSRP Summary of this function goes here
%   Detailed explanation goes here

serialNumberRx = '30C51BC';
serialNumberTx = '30C5426';

devicesToCheck = [];
if strcmp(mode, 'loopback')
    serialNumberRx = serialNumberTx;
    if samplesPerFrame > 1e5
        warning('samplePerFrame > 1e5 leads to overrun: the board drops samples.')
    end
    devicesToCheck = [serialNumberRx];
elseif strcmp(mode, 'oneBoardTx')
    devicesToCheck = [serialNumberTx];
elseif strcmp(mode, 'oneBoardRx')
    devicesToCheck = [serialNumberRx];
elseif strcmp(mode, 'twoBoardsRxTx')
    devicesToCheck = [serialNumberRx; serialNumberTx];
end

checkDevicesConnection(devicesToCheck);

boardPlatform = 'B200';
fc = 2.4e9; % carrier frequency
LOO = 1e3; % local oscillator offset
clockRateTx = 5e6; % main clock (Sampling rate of the digital signal sent to ADC) (5e6 to 56e6)
interpolationTx = 10; % Interpolation factor from host signal to USRP signal (e.g Fs = 0.5 Mhz => ClockRate = 5 Mhz)
clockRateRx = clockRateTx; % on the same radio we need to use the same
decimationRx = interpolationTx; % to get to 1MHz
clockInputSource = 'Internal';
outputDataTypeUSRP = 'double'; % Difference between transport and output data type ?
gainTx = 60; %60 loopback 30dB attenuator; % 89 over the air, no attenuator
gainRx = 40;

burstMode = false;

Ts = interpolationTx/clockRateTx; % [s] sampling time
Fs = 1/Ts;

noFrames = length(dataTx)/samplesPerFrame;
if mod(noFrames, 1) ~= 0
    error('Number of frames must be integer %s.', num2str(noFrames));
end

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

if strcmp(mode, 'loopback') || strcmp(mode, 'oneBoardRx') ||  strcmp(mode, 'twoBoardsRxTx')
    
    dataRx = zeros(1, noFrames*samplesPerFrame);
    
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

switch mode
    case {'twoBoardsRxTx','loopback'}
        l = 1;
        while l <= noFrames
            [dataRxUSRP, len, ~] = rx();
            
            if len > 0
                dataRx(1+(l-1)*samplesPerFrame:l*samplesPerFrame) = dataRxUSRP.';
                tx(dataTx(1+(l-1)*samplesPerFrame:l*samplesPerFrame).');
                l = l+1;
            end
        end
        
    case 'oneBoardRx'
        l = 1;
        while l <= noFrames
            [dataRxUSRP, len, ~] = rx();
            
            if len > 0
                dataRx(1+(l-1)*samplesPerFrame:l*samplesPerFrame) = dataRxUSRP.';
                l = l+1;
            end
        end
    case 'oneBoardTx'
        for l = 1:noFrames
            tx(dataTx(1+(l-1)*samplesPerFrame:l*samplesPerFrame).');
            dataRx = [];
        end
    otherwise
        error('Unknown mode of operation.')
end

if strcmp(mode, 'loopback') || strcmp(mode, 'oneBoardTx') ||  strcmp(mode, 'twoBoardsRxTx')
    release(tx);
end

if strcmp(mode, 'loopback') || strcmp(mode, 'oneBoardRx') ||  strcmp(mode, 'twoBoardsRxTx')
    release(rx);
end   

end

