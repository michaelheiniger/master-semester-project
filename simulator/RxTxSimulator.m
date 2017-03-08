classdef RxTxSimulator < TxRxInterface
    %RXTXSIMULATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = private, GetAccess = public)
        Fs % Hz
        SNR % dB
        awgnEnabled % 0 or 1
        doppler % Hz, amount of Doppler shift
        clockOffsetEnabled % 0 or 1
        USFOffsetSimulation % integer, simulate offset between Rx/Tx clocks
        clockSampleOffset % Between 1 and USFOffsetSimulation-1
        randomCutEnabled % 0 or 1
        startPos % integer
    end
    
    properties (Dependent)
        Ts
    end
    
    methods (Access = public)
        function obj = RxTxSimulator(lengthPulse)
            obj.Fs = 0.5e6;
            
            % AWGN
            obj.SNR = 30;
            obj.awgnEnabled = 1;
            
            % Doppler
            obj.doppler = 234;
            
            % Clock offset
            obj.clockOffsetEnabled = 1;
            obj.USFOffsetSimulation = 4;
            obj.clockSampleOffset = 2;
            
            % Random cut
            obj.randomCutEnabled = 1;
            obj.startPos = lengthPulse+123;
        end
        
        function dataRx = transmitReceive(obj, dataTx, ~)
            
            dataRx = dataTx;
            
            %Add Complex White Gaussian Noise
            if obj.awgnEnabled == 1
                dataRx = awgn(complex(dataRx), obj.SNR, 'measured');
            end
            
            % Add constant Doppler shift
            if obj.doppler ~= 0
                fprintf('Doppler added: %f\n', obj.doppler);
                t = (0:(length(dataRx)-1))*obj.Ts;
                dataRx = dataRx.*exp(1j*2*pi*obj.doppler*t);
            else
                disp('No Doppler');
            end
            
            % Fixed clock offset simulation
            if obj.clockOffsetEnabled == 1
                disp('Clock offset added');
                dataRx = resample(dataRx, obj.USFOffsetSimulation, 1);
                dataRx = downsample(dataRx, obj.USFOffsetSimulation, obj.clockSampleOffset);
            else
                disp('No clock offset added');
            end
            
            % Cut signal at random position
            if obj.randomCutEnabled == 1
                disp('Received signal cut');
                dataRx = dataRx(obj.startPos:end);
            else
                disp('Received signal NOT cut');
            end
        end
        
        function releaseTransmitter(obj)
            % Nothing to do
        end
        
        function releaseReceiver(obj)
            % Nothing to do
        end
    end
    
    methods
        function Ts = get.Ts(obj)
            Ts = 1/obj.Fs;
        end
    end
    
end

