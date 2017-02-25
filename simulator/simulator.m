%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulator: ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulator parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data
load('CAcodes.mat');
ca = satCAcodes(1,:); % Data
N = 10; % Number of repetition of CA code (N > 1)
symbols = repmat(ca,1,N);
numFullCAToKeep = 2; % (> 1) Number of CA to keep (first complete one included)

% Pulse shaping
Fd = 1;         % symbol rate [symbols/s]
USF = 10;       % upsampling factor
Fs = USF*Fd;    % sampling frequency
Ts = 1/Fs;

fprintf('Fs: %f, Ts: %f\n', Fs, Ts);

% We need to truncate the infinite impulse response of the RRC filter
% to approximate it by a FIR filter
delay = 10;  % in number of symbol periods
			 % ( = number of lobes of the time response at each side of the main one that we keep).
beta = 0.5;  % roll-off factor
pulse = rcosdesign(beta, 2*delay, Fs, 'sqrt'); % already normalized to 1
% fvtool(pulse, 'Analysis', 'impulse')   % Visualize the filter

% pulse1 = [pulse, zeros(1,Fd*USF)];
% pulse2 = [zeros(1,Fd*USF), pulse];
% subplot(1,2,1),plot(pulse1,'.');
% hold on;
% plot(pulse2,'.');
% subplot(1,2,2),plot(pulse1+pulse2);

%
caUpsampled = upsample(ca,USF);

% AWGN
SNR = 5; % Signal-To-Noise ratio in dB
awgnEnabled = 1;

% Doppler
dopplerCorrection = 1;
doppler = 0.000132;%0.0001 % Amount of Doppler for simulation
maxDoppler = 0.001; % Absolute value defining the range in which to estimate the Doppler
% dopplerStep = 0.1; % Size of the steps in the range of Doppler
dopplerStep = 0.1*2*pi/((length(caUpsampled)-1)*Ts);
fprintf('DopplerStep: %f\n',dopplerStep);
fcorr = 0.0001;

% Clock offset
clockOffsetEnabled = 1;
clockOffsetCorrection = 1;
USFOffsetSimulation = 4; % USF used to simulate clock offset
clockSampleOffset = 1; % Offset (in #samples) of the clock
USFClockOffsetCorrection = 10; % USF used to correct for the clock offset

% Random cut
randomCutEnabled = 1;
startPos = randi([1,length(caUpsampled)/2]);

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmitter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate symbol-by-symbol pulse train samples
txSignal = symbolsToSamples(symbols, pulse, Fs);

rxSignal = txSignal;
rxSignalImpaired = rxSignal;
% length(tx_signal) % length(pulse)+#symbols*USF-1 

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Impairement simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rxSignalImpaired = [zeros(1,1000), rxSignalImpaired, zeros(1,1000)];    

%Add Complex White Gaussian Noise
% complex() makes sure that the signal is considered to be complex so that
% awgn() adds complex noise.
if awgnEnabled == 1
    rxSignalImpaired = awgn(complex(rxSignalImpaired), SNR, 'measured');
end

% Add constant Doppler shift
t = (0:(length(rxSignalImpaired)-1))*Ts;
if doppler ~= 0
    fprintf('Doppler added: %f\n', doppler);
    rxSignalImpaired = rxSignalImpaired.*exp(1j*2*pi*doppler*t);
else
   disp('No Doppler'); 
end

% Fixed clock offset simulation
if clockOffsetEnabled == 1
    disp('Clock offset added'); 
    rxSignalImpaired = resample(rxSignalImpaired, USFOffsetSimulation, 1);
%     rxSignalImpaired = interp(rxSignalImpaired, USFOffsetSimulation);
    rxSignalImpaired = downsample(rxSignalImpaired, USFOffsetSimulation, clockSampleOffset);
else
    disp('No clock offset'); 
end

% Cut signal at random position
if randomCutEnabled == 1
    disp('Received signal cut'); 
    rxSignalImpaired = rxSignalImpaired(startPos:end);
else
    disp('Received signal NOT cut'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matched-Filter
matchedFilter = conj(fliplr(pulse));
matchedFilterOutput = conv(rxSignalImpaired, matchedFilter);

% Doppler shift coarse estimation
[tauEstimate, coarseDopplerEstimate] = dopplerCoarseEstimate(matchedFilterOutput, Ts, caUpsampled, maxDoppler, dopplerStep);
fprintf('Tau: %d, Doppler coarse: %f, doppler range: [%f,%f], step: %f\n', tauEstimate,coarseDopplerEstimate, -maxDoppler/2, maxDoppler/2, dopplerStep)

% Doppler shift fine estimation
if dopplerCorrection == 1
    disp('Doppler correction');
    fineDopplerEstimate = dopplerFineEstimate(matchedFilterOutput(tauEstimate:end), caUpsampled, Ts, coarseDopplerEstimate, fcorr);
    fprintf('Fine doppler estimate: %f\n', fineDopplerEstimate);
    
    % Doppler shift removal
    t = (0:(length(matchedFilterOutput)-1))*Ts;
    matchedFilterOutput = matchedFilterOutput .* exp(-1j*2*pi*fineDopplerEstimate*t);
    
    % Intial phase estimation
    initialPhaseEstimate = angle(sum(matchedFilterOutput(tauEstimate:tauEstimate+length(caUpsampled)-1).*conj(caUpsampled)));
    
    % Initial phase removal
    matchedFilterOutput = matchedFilterOutput .* exp(-1j*initialPhaseEstimate);
    
%     [R3, ~] = xcorr(matchedFilterOutput, caUpsampled);
%     R3 = R3(length(matchedFilterOutput):end);
%     [~, tauEstimateTmp] = max(abs(R3))
%     plot(R3);
    
end

% Correct for clock offset by upsampling and correlating
if clockOffsetCorrection == 1
    disp('Clock offset correction');
    caUpsampled = upsample(caUpsampled, USFClockOffsetCorrection);
    matchedFilterOutput = resample(matchedFilterOutput, USFClockOffsetCorrection, 1);
%     matchedFilterOutput = interp(matchedFilterOutput, USFClockOffsetCorrection);
    [R2, ~] = xcorr(matchedFilterOutput, caUpsampled);
    R2 = R2(length(matchedFilterOutput):end);
    [~, tauEstimate] = max(abs(R2));
    fprintf('New Tau: %d\n', tauEstimate);
%     plot(R2);
else
    USFClockOffsetCorrection = 1;
end

% Sample the matched filter to get the symbols estimates
% symbolsEstimates = matchedFilterOutput(tauEstimate:USF*USFClockOffsetCorrection:end-length(pulse)*USFClockOffsetCorrection+1);
symbolsEstimates = matchedFilterOutput(tauEstimate:USF*USFClockOffsetCorrection:tauEstimate+length(caUpsampled)*numFullCAToKeep-1);

% Plot the constellation
scatter(real(symbolsEstimates), imag(symbolsEstimates),'.');
hold on;
scatter([-1 1],[0 0], 'kx');
axis([-1.5,1.5,-1.5,1.5]);
grid on;

map = pskMap(2);
% size(symbols(1:length(ca)*numFullCAToKeep))
% size(symbolsEstimates)
symbolsDec = demodulator(symbols(1:length(ca)*numFullCAToKeep),map);
symbolsEstimatesDec = demodulator(symbolsEstimates,map);
BER = sum(symbolsDec ~= symbolsEstimatesDec)/length(symbolsDec)

