function [vector] = changeVectorBackToRowIfNeeded(vector, wasColumnVector)
%CHANGEVECTORBACKTOROWIFNEEDED change back the vector to a row-vector if
% needed.

numRows = size(vector, 1);
numCols = size(vector, 2);

if numRows > 1 && numCols > 1
    error(['The input is NOT a vector: ', num2str(size(vector,1)), 'x', num2str(size(vector,2))]);
end

if not(wasColumnVector)
    vector = transpose(vector);
end

end

