function [tauEstimate, dopplerEstimate] = dopplerCoarseEstimate(signal, Ts, referenceSignal, maxDoppler, dopplerStep)
%DOPPLERCOARSEESTIMATE Summary of this function goes here
%   Detailed explanation goes here

if length(signal) < length(referenceSignal)
    error('Signal must at least as long as the reference signal.');
elseif length(signal) >= 2*length(referenceSignal)-1
    signalTruncated = signal(1:2*length(referenceSignal)-1);
    disp('Doppler coarse case 2');
else
    signalTruncated = signal;
    disp('Doppler coarse case 3');
end

% Time axis used to correct for Doppler shift
t = (0 : length(signalTruncated)-1)*Ts;

% Range of Doppler shifts to try
dopplerRange = -maxDoppler/2:dopplerStep:maxDoppler/2;

% Save the value of the max correlation
maxCorrelation = 0;

tauEstimate = 1;

% Exhaustive search of the best Doppler shift in the given range
for k = dopplerRange
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correct for Doppler: Remove Doppler shift 
    % This needs to be done before computing the cross-correlation
    % in order to optimize the results since the reference signal and 
    % the received signal have different frequency due to Doppler shift
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    correctedSignal = exp(-1j*2*pi*k*t) .* signalTruncated;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cross-correlation of reference signal with received signal to find
    % tau and synchronize on the sample level
    % max() must be used since the CA code is multiplied by -1 or +1
    % and the reference sequence corresponds to +1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [R, ~] = xcorr(correctedSignal, referenceSignal);
    
    % Tau cannot be < 0 since the max must be when overlap is complete
    R = R(length(correctedSignal):end);
    R = R(1:length(referenceSignal));
    
    [correlation, tau] = max(abs(R));
    
    if correlation > maxCorrelation
        dopplerEstimate = k;
        tauEstimate = tau;
        maxCorrelation = correlation;
%         plot(R)
%         pause
    end    
end

% Take modulo here to get the start of the first CA
tauEstimate = mod(tauEstimate, length(referenceSignal));
% display(['Tau estim: ' num2str(best_tau)]);

% correctedSignal = exp(-1j*2*pi*dopplerEstimate*t) .* signalTruncated;
% initialPhaseEstimate = angle(sum(correctedSignal(tauEstimate:tauEstimate+length(referenceSignal)-1).*conj(referenceSignal)));
% initialPhaseEstimate = angle(sum(correctedSignal(1:length(referenceSignal)).*conj(referenceSignal)));

end


