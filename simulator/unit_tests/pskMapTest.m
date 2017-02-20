function [ result ] = pskMapTest()
%PSKMAPTEST Summary of this function goes here
%   Detailed explanation goes here

t1 = isequal(pskMap(2), [1.0000 + 0.0000i  -1.0000 + 0.0000i]);

t2 = isequal(pskMap(4), [1.0000 + 0.0000i   0.0000 + 1.0000i  -1.0000 + 0.0000i  -0.0000 - 1.0000i]);

t3 = isequal(pskMap(8), [1.0000 + 0.0000i   0.7071 + 0.7071i   0.0000 + 1.0000i  -0.7071 + 0.7071i -1.0000 + 0.0000i  -0.7071 - 0.7071i  -0.0000 - 1.0000i   0.7071 - 0.7071i]);

result = mod(t1+t2+t3,3) == 0;

end

