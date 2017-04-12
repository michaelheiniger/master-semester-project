function [innerProduct] = doInnerProduct(v1, v2)
%INNERPRODUCT Does the inner product of two vectors: v1 * conj(v2)

if not(isvector(v1)) || not(isvector(v2))
   error('Inputs must be vectors (column [n 1] or row [1 n]).'); 
end

% Make v2 a column-vector
v2 = v2(:);

[~, col] = size(v1);
if col == 1
    v1 = transpose(v1);
end

innerProduct = v1 * conj(v2);

end

