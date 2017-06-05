function [vector, wasColumnVector] = turnIntoColumnVector(vector)
%TURNINTOCOLUMNVECTOR returns the vector as a column-vector and
% wasColumnVector == 1 if it was already a column-vector and 0 if it was a
% row-vector

numRows = size(vector, 1);
numCols = size(vector, 2);

if numRows > 1 && numCols > 1
    error(['The input is NOT a vector: ', num2str(size(vector,1)), 'x', num2str(size(vector,2))]);
end

isRowVector = numCols > 1;

wasColumnVector = 1;
if isRowVector
    wasColumnVector = 0;
    vector = transpose(vector);
end

end

