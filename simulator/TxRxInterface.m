classdef TxRxInterface < handle
    %TXRXINTERFACE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Abstract)
        dataRx = transmitReceive(obj, dataTx, noFramesRx);
        
        releaseTransmitter(obj);
        
        releaseReceiver(obj);
    end
    
end

