function [] = checkDevicesConnection(serialNumbers)
%CHECKDEVICESCONNECTION Check if the devices are connected to the host and
%ready to be used.
% [] = checkDevicesConnection(serialNumbers): serialNumbers is a matrix
% containing one serialNumber on every row. Duplicate rows are removed.

serialNumbersToCheck = unique(serialNumbers, 'rows');

for k = 1:size(serialNumbersToCheck,1)
    usrp = findsdru(serialNumbersToCheck(k,:));
    
    if not(strcmp(usrp.Status,'Success'))
        error('Connection error with the device.');
    end
end
end

