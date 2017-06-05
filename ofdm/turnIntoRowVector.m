function [vector, wasRowVector] = turnIntoRowVector(vector)
%TURNINTOROWVECTOR returns the vector as a row-vector and
% wasRowVector == 1 if it was already a row-vector and 0 if it was a
% coulmn-vector

numRows = size(vector, 1);
numCols = size(vector, 2);

if numRows > 1 && numCols > 1
    error(['The input is NOT a vector: ', num2str(size(vector,1)), 'x', num2str(size(vector,2))]);
end

isColumnVector = numRows > 1;

wasRowVector = 1;
if isColumnVector
    wasRowVector = 0;
    vector = transpose(vector);
end

end

