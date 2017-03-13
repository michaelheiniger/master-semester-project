classdef RxTxUSRP < TxRxInterface
    %RXTXUSRP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private, GetAccess = public)
        
        boardPlatform
        serialNumberTx
        serialNumberRx
        rxEnabled
        txEnabled
        tx
        rx
        
        fc % carrier frequency
        LOO % local oscillator offset
        clockRate % main clock (Sampling rate of the digital signal sent to ADC) (5e6 to 56e6)
        interpolationTx % Interpolation factor from host signal to USRP signal (e.g Fs = 0.5 Mhz => ClockRate = 5 Mhz)
        decimationRx
        clockInputSource
        outputDataTypeUSRP % Difference between transport and output data type ?
        gainTx
        gainRx
        samplesPerFrame
        burstMode
        noFramesInBurst
    end
    
    properties (Dependent, GetAccess = public)
        Ts
    end
    
    methods (Access = public)
        function obj = RxTxUSRP(boardPlatform, serialNumberRx, serialNumberTx, rxEnabled, txEnabled, burstMode, samplesPerFrame)
            
            if nargin == 0
                error('Missing arguments');
            end
            
            obj.boardPlatform = boardPlatform;
            obj.serialNumberRx = serialNumberRx;
            obj.serialNumberTx = serialNumberTx;
            obj.rxEnabled = rxEnabled;
            obj.txEnabled = txEnabled;
            obj.fc = 2.4e9;
            obj.LOO = 100e3;
            obj.clockRate = 5e6;
            obj.interpolationTx = 10;
            obj.decimationRx = obj.interpolationTx;
            obj.clockInputSource = 'Internal';
            obj.outputDataTypeUSRP = 'double';
            obj.gainTx = 60; %60 loopback 30dB attenuator; % 89 over the air, no attenuator
            obj.gainRx = 40;
            obj.samplesPerFrame = samplesPerFrame;
            obj.burstMode = burstMode;
            obj.noFramesInBurst = 154;
            
            if rxEnabled
                obj.rx = getReceiver(obj);
            end
            if txEnabled
                obj.tx = getTransmitter(obj);
            end
            
        end
        
        function dataRx = transmitReceive(obj, dataTx, noFramesRx)
            
            if mod(length(dataTx),obj.samplesPerFrame) ~= 0
                error('The number of samples to transmit is not a multiple of the frame size.');
            end
            
            noFramesTx = ceil(length(dataTx)/obj.samplesPerFrame);
            
            obj.checkDevicesConnection();
            
            dataRx = zeros(1, noFramesRx*obj.samplesPerFrame);
            
            k = 1;
            l = 1;
            while k <= noFramesRx || l <= noFramesTx
                
                if k <= noFramesRx && obj.rxEnabled
                    % Receive samples
                    [frameRx, len, ~] = obj.rx();
                    
                    if len > 0
                        % Save received samples
                        dataRx(1+(k-1)*obj.samplesPerFrame:k*obj.samplesPerFrame) = frameRx.';
                        k = k+1;
                    end
                end
                
                if l <= noFramesTx && obj.txEnabled && not(isempty(length(dataTx)))
                    % Transmit samples
                    obj.tx(dataTx(1+(l-1)*obj.samplesPerFrame:l*obj.samplesPerFrame).');
                    l = l+1;
                end
            end
        end
        
        function releaseTransmitter(obj)
            if obj.txEnabled
                release(obj.tx);
            end
        end
        
        function releaseReceiver(obj)
            if obj.rxEnabled
                release(obj.rx);
            end
        end
        
    end
    methods
        function Ts = get.Ts(obj)
            Ts = obj.interpolationTx/obj.clockRate;
        end
    end
    
    methods (Access = private)
         
        function tx = getTransmitter(obj)
            
            tx = comm.SDRuTransmitter(...
                'Platform',obj.boardPlatform, ...
                'CenterFrequency', obj.fc, ...
                'LocalOscillatorOffset', obj.LOO,...
                'InterpolationFactor', obj.interpolationTx, ...
                'MasterClockRate', obj.clockRate,...
                'SerialNum', obj.serialNumberTx,...
                'Gain', obj.gainTx,...
                'ClockSource', obj.clockInputSource,...
                'UnderrunOutputPort', true);
            
            if obj.burstMode
                tx.EnableBurstMode = true;
                tx.NumFramesInBurst = obj.noFramesInBurst;
            end
            
        end
        
        function rx = getReceiver(obj)
            
            rx = comm.SDRuReceiver(...
                'Platform', obj.boardPlatform, ...
                'CenterFrequency', obj.fc, ...
                'LocalOscillatorOffset', obj.LOO,...
                'MasterClockRate', obj.clockRate,...
                'SerialNum', obj.serialNumberRx,...
                'DecimationFactor', obj.decimationRx,...
                'Gain', obj.gainRx,...
                'SamplesPerFrame', obj.samplesPerFrame,...
                'ClockSource', obj.clockInputSource,...
                'OutputDataType', obj.outputDataTypeUSRP,...
                'OverrunOutputPort', true);
            
            if obj.burstMode
                rx.EnableBurstMode = true;
                rx.NumFramesInBurst = obj.noFramesInBurst;
            end
        end
        
        function [] = checkDevicesConnection(obj)
                        
            if strcmp(obj.serialNumberRx, obj.serialNumberTx)
                if obj.rxEnabled
                    usrpRx = findsdru(obj.serialNumberRx);
                    usrpTx = usrpRx;
                elseif obj.txEnabled
                   usrpTx = findsdru(obj.serialNumberTx); 
                   usrpRx = usrpTx;
                end
            else
                if obj.rxEnabled
                    usrpRx = findsdru(obj.serialNumberRx);
                end
                if obj.txEnabled
                    usrpTx = findsdru(obj.serialNumberTx); 
                end
            end
            
            if not(strcmp(usrpRx.Status,'Success')) || not(strcmp(usrpTx.Status,'Success'))
                error('Connection error with the device(s).');
            end
        end
    end
end

