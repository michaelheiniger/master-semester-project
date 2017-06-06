function [timingEst, frequencyOffsetEst, initialPhaseEst] = timingAndFrequencyOffsetEstimation2(systemConfig, receiverConfig, signal, ca)
%TIMINGANDFREQUENCYOFFSETESTIMATION2 Find the beginning of the frame, the
%frequency offset and intial phase.
% Returned values are:
% - timingEst: estimate of the index of the first sample of the frame
% - frequencyOffsetEst: time-dependent frequency offset due different to
% carrier frequency at transmitter and receiver
% - intialPhaseEst: estimate intial phase offset (i.e. does not depend on time)
% NOTE: "Doppler" is used here as synonym for "frequency offset" to be
% consistent with the GPS part of the SDR course since this is where this 
% technique comes from

sc = systemConfig;
rc = receiverConfig;

Ts = 1/sc.Fs;

maxDoppler = 1000;

% Doppler step computation
dopplerStep = 0.01/((length(ca)-1)*Ts);
fprintf('DopplerStep: %f\n', dopplerStep);

%%%%%%%%%%%%%%%%%%%%%%%%% START COARSE EST %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time axis used to correct for Doppler shift
t = (0 : length(signal)-1)*Ts;

% Range of Doppler shifts to try
dopplerRange = -maxDoppler/2:dopplerStep:maxDoppler/2;

% Save the value of the max correlation
maxCorrelation = 0;

tauEstimate = 1;

% Exhaustive search of the best Doppler shift in the given range
for k = dopplerRange
    
    correctedSignal = exp(-1j*2*pi*k*t) .* signal;
        
    [R, ~] = xcorr(correctedSignal, ca);
    
    % Tau cannot be < 0 since the max must be when overlap is complete
    R = R(length(correctedSignal):end);
    
    [correlation, tau] = max(abs(R));
    
    if correlation > maxCorrelation
        coarseDopplerEstimate = k;
        tauEstimate = tau;
        maxCorrelation = correlation;
    end    
end

disp(['Tau est: ', num2str(tauEstimate)]);
disp(['Coarse Doppler est: ', num2str(coarseDopplerEstimate)]);

% Tau fine estimation
ISILength = 30;
results = zeros(ISILength,1);
for i = 0:ISILength-1
    tentativeTau = tauEstimate-i;
    caRx = signal(tentativeTau:tentativeTau+length(ca)-1).';
    h = channel_estimation_causal(caRx,ca);
%     if tentativeTau == 1001
        figure;
        stem(real(h),'.');
        title(['Tentative tau: ' num2str(tentativeTau)]);
%     end
    results(i+1) = sum(abs(real(h(1:ISILength))));
%     results(i+1) = sum(abs(real(h(1:16))));
end

figure;
plot(tauEstimate-(0:ISILength-1),results);
% plot(results);
title('abs sum of first 16 taps of h vs tentative tau');



%%%%%%%%%%%%%%%%%%%%%%%%% Upsampled version %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% signalUp = resample(signal, rc.USF,1);
% caUp = upsample(ca, rc.USF);
% 
% % Time axis used to correct for Doppler shift
% t = (0 : length(signalUp)-1)*Ts;
% 
% % Range of Doppler shifts to try
% dopplerRange = -maxDoppler/2:dopplerStep:maxDoppler/2;
% 
% % Save the value of the max correlation
% maxCorrelationUp = 0;
% 
% tauEstimateUp = 1;
% 
% % Exhaustive search of the best Doppler shift in the given range
% for k = dopplerRange
%     
%     correctedSignal = exp(-1j*2*pi*k*t) .* signalUp;
%         
%     [R, ~] = xcorr(correctedSignal, caUp);
%     
%     % Tau cannot be < 0 since the max must be when overlap is complete
%     R = R(length(correctedSignal):end);
% %     R = R(1:length(ca));
%     
%     [correlation, tau] = max(abs(R));
%     
%     if correlation > maxCorrelationUp
%         coarseDopplerEstimateUp = k;
%         tauEstimateUp = tau;
%         maxCorrelationUp = correlation;
%     end    
% end
% 
% disp(['Tau est (upsampled): ', num2str(tauEstimateUp)]);
% disp(['Coarse Doppler est (upsampled): ', num2str(coarseDopplerEstimateUp)]);

%%%%%%%%%%%%%%%%%%%%%%%%% END COARSE EST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% START FINE EST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

caRx = signal(tauEstimate:tauEstimate+length(ca)-1);

dopplerRange = [coarseDopplerEstimate-dopplerStep, coarseDopplerEstimate, coarseDopplerEstimate+dopplerStep]; 
innerProducts = zeros(length(dopplerRange), 1);

t = (0:(length(ca)-1))*Ts;
for k = 1:length(dopplerRange)
    caRxCorrected = caRx .* exp(-1j*2*pi*dopplerRange(k)*t);
    innerProducts(k) = sum(ca .* conj(caRxCorrected));
end

d1 = coarseDopplerEstimate-dopplerStep;
d2 = coarseDopplerEstimate;
d3 = coarseDopplerEstimate+dopplerStep;

D = [ d1.^2, d1, 1 ; 
   d2.^2, d2, 1; 
   d3.^2, d3, 1 ];

f = [ innerProducts(1); innerProducts(2) ;innerProducts(3) ]; 
result = D\abs(f);
a = result(1);
b = result(2);

fineDopplerEstimate = -b/(2*a);

disp(['Fine Doppler est: ', num2str(fineDopplerEstimate)])

%%%%%%%%%%%%%%%%%%%%%%%%% END FINE EST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timingEst = tauEstimate;
frequencyOffsetEst = fineDopplerEstimate;

% Compute initial phase estimate
t = (0:(length(ca)-1))*Ts;
caRxCorrected = caRx .* exp(-1j*2*pi*fineDopplerEstimate*t);
initialPhaseEst = angle(doInnerProduct(caRxCorrected, ca));

end

