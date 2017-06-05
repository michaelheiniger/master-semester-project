function [isReceiver] = isInstanceReceiver(mode)
%ISINSTANCERECEIVER determines if the instance is a receiver as a
%function of the mode of operation

if strcmp(mode, 'simulation') || strcmp(mode, 'twoBoardsRxTx') || strcmp(mode, 'oneBoardRx') || strcmp(mode, 'loopback')
    isReceiver = 1;
else
    isReceiver = 0;
end

end

