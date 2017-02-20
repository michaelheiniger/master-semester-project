function [doppler_estim, tau_estim] = dopplerCoarseEstimate(data, fc, Ts, USF, maxDoppler, dopplerStep)
%DOPPLERCOARSEESTIMATE Summary of this function goes here
%   Detailed explanation goes here

load('CAcodes.mat');
ca = satCAcodes(1,:);
ca_us = upsample(ca,USF);

dataTruncated = data; %(1:2*length(ca)-1);

% Save the tau corresponding to the max correlation
best_tau = 0;

% Time axis used to correct for Doppler shift
t = (0 : length(dataTruncated)-1)*Ts;

% Range of Doppler shifts to try
dopplerShift = -maxDoppler/2:dopplerStep:maxDoppler/2;

% Save the value of the max correlation
bestCorrelation = 0;

% Save the value of the Doppler shift that gives the max correlation
bestDopplerShift = 0;

% Exhaustive search of the best Doppler shift in the given range
for k = dopplerShift
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correct for Doppler: Remove Doppler shift 
    % This needs to be done before computing the cross-correlation
    % in order to optimize the results since the reference signal and 
    % the received signal have different frequency due to Doppler shift
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    shifted_data = exp(-1j*2*pi*fc*k*t).*dataTruncated;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cross-correlation of reference signal with received signal to find
    % tau and synchronize on the sample level
    % max() must be used since the CA code is multiplied by -1 or +1
    % and the reference sequence corresponds to +1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Data part 1
    [R, ~] = xcorr(shifted_data, ca_us);
    
    % Tau cannot be < 0 since the max must be when overlap is complete
    R = R(length(shifted_data):end);
    
    % Tau cannot be > Ta*fs where Ta is the time of a full CA code
    R = R(1:length(ca_us));
    [max_corr_data, tau] = max(abs(R));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If the resulting value is higher than the maximum correlation so far,
    % the best value is saved and the tau and Doppler shifts estimates are
    % updated as well.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if max_corr_data > bestCorrelation
        bestDopplerShift = k;
        best_tau = tau;
        bestCorrelation = max_corr_data;
    end    
end

doppler_estim = bestDopplerShift;
% Take modulo here to get the start of the first CA
tau_estim = mod(best_tau, length(ca_us));

% display(['Tau estim: ' num2str(best_tau)]);

end


