function [bits] = getBitsToSend(source, numBits, dataSourceConfig)
%GETBITSTOSEND Summary of this function goes here
%   Detailed explanation goes here

dsc = dataSourceConfig;

switch source
    
    case 'random'
        bits = randi([0,1], numBits, 1);
        
    case 'textFile'
        bitsRead = textFileToBits(dsc.filePath);
        
        bitsReadLength = length(bitsRead);
        if bitsReadLength < numBits
           numRep = ceil(numBits/bitsReadLength);
           bitsRead = repmat(bitsRead, numRep, 1);
        end
        
        % NOTE: numBits may not be a multiple of 8,
        % Text recovery must truncate the last k bits in this case
        % where k is in {1,...,7}
        bits = bitsRead(1:numBits);
        
    otherwise
        error('getBitsToSend: source of bits unknown')

end

