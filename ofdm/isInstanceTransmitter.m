function [isTransmitter] = isInstanceTransmitter(mode)
%ISINSTANCETRANSMITTER determines if the instance is a transmitter as a
%function of the mode of operation

if strcmp(mode, 'simulation') || strcmp(mode, 'oneBoardTx') || strcmp(mode, 'loopback') || strcmp(mode, 'twoBoardsRxTx')
    isTransmitter = 1;
else
    isTransmitter = 0;
end

end

