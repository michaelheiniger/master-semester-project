function [] = unitTestResult( testNumber, testResult )
%UNITTESTERROR Summary of this function goes here
%   Detailed explanation goes here

if testResult == 1
    sprintf('Test %d succeeded',testNumber)
else
    sprintf('Test %d failed',testNumber)
end


end

