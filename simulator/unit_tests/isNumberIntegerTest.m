function [ result ] = isNumberIntegerTest()
%ISINTEGERTEST Summary of this function goes here
%   Detailed explanation goes here

t1 = not(isNumberInteger(2.5));

t2 = isNumberInteger(3);

t3 = not(isNumberInteger(pi));

result = mod(t1+t2+t3,3) == 0;

end

