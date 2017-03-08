function [ signalCorrected, tauEst ] = dopplerEstimationAndCorrection(signal, referenceSignal, Ts, maxDoppler, dopplerCorrection)
%DOPPLERESTIMATIONANDCORRECTION Summary of this function goes here
%   Detailed explanation goes here

% Doppler step computation
dopplerStep = 0.1/((length(referenceSignal)-1)*Ts);
fprintf('DopplerStep: %f\n', dopplerStep);

% Doppler shift coarse estimation
[tauEst, dopplerCoarseEst] = dopplerCoarseEstimate(signal, Ts, referenceSignal, maxDoppler, dopplerStep);
fprintf('Tau: %d, Coarse doppler estimate: %f, doppler range: [%f,%f], step: %f\n', tauEst, dopplerCoarseEst, -maxDoppler/2, maxDoppler/2, dopplerStep)

% Doppler shift fine estimation
dopplerFineEst = dopplerFineEstimate(signal(tauEst:end), referenceSignal, Ts, dopplerCoarseEst, dopplerStep);
fprintf('Fine doppler estimate: %f\n', dopplerFineEst);

% Remove irrelevant part of the signal
signalCorrected = signal(tauEst:end);

% Doppler shift correction
if dopplerCorrection == 1
    disp('Doppler correction');
    
    % Doppler shift removal
    t = (0:(length(signalCorrected)-1))*Ts;
    signalCorrected = signalCorrected .* exp(-1j*2*pi*dopplerFineEst*t);
    
    % Intial phase estimation
    initialPhaseEst = initialPhaseEstimate(signalCorrected, referenceSignal);
    fprintf('Initial phase estimate: %f\n', initialPhaseEst);
    
    % Initial phase removal
    signalCorrected = signalCorrected .* exp(-1j*initialPhaseEst);
end

end
