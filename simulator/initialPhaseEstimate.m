function [ initialPhaseEst ] = initialPhaseEstimate(signal, referenceSignal )
%INITIALPHASEESTIMATE Summary of this function goes here
%   Detailed explanation goes here

    initialPhaseEst = angle(sum(signal(1:length(referenceSignal)).*conj(referenceSignal)));


end

