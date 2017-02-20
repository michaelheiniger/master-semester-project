function [ result ] = modulatorTest()
%MODULATORTEST Summary of this function goes here
%   Detailed explanation goes here

% TO WRITE !

x = [0 1 2 3 2 2 0 1 3];
t1 = isequal(modulator(x, map(4)), [1.0000 + 0.0000i   0.0000 + 1.0000i  -1.0000 + 0.0000i  -0.0000 - 1.0000i  -1.0000 + 0.0000i  -1.0000 + 0.0000i   1.0000 + 0.0000i   0.0000 + 1.0000i -0.0000 - 1.0000i]);

xMat = [ 0 1 2 3; 1 2 3 4];
t2 = modulator(xMat, map(4));

t3 = 

result = mod(t1+t2+t3,3) == 0;

end

