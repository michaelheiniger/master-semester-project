function [interpolationFactor] = adjustUsrpInterpolationFactor(Fs)
%ADJUSTUSRPINTERPOLATIONFACTOR Summary of this function goes here
%   Detailed explanation goes here

if Fs < 5e6
    interpolationFactor = 10;
else
    interpolationFactor = 1;
end

end

