clear;
clc;

% Signal-to-Noise ratio
SNR = 10; % in dB

% Number of repetition of CA code
N = 4;

% Get data
load('CAcodes.mat');
ca = satCAcodes(1,:);
preamble = repmat(ca,1,N);
% preamble = preamble(1:4);

% Translate preamble -1 -> 1 and 1 -> 0
preambleBits = (preamble+1)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmitter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bitsPerFrame = 2000;
payload = []; %randi([0,1],1,bitsPerFrame); % Information bits 
data = [preambleBits payload];

% BSPK Modulation
map = pskMap(2);
symbols = modulator(data, map);

% Pulse shaping
Fd = 1;      % symbol rate [symbols/s]
USF = 10;     % upsampling factor
Fs = USF*Fd % sampling frequency
Ts = 1/Fs

% We need to truncate the infinite impulse response of the RRC filter 
% to approximate it by a FIR filter
delay = 10;  % in number of symbol periods 
			 % ( = number of lobes of the time response at each side of the main one that we keep).

beta = 0.5; % roll-off factor 

pulse = rcosdesign(beta, 2*delay, Fs, 'sqrt'); % already normalized to 1
% fvtool(pulse, 'Analysis', 'impulse')   % Visualize the filter

% Generate symbol-by-symbol pulse train samples
tx_signal = symbolsToSamples(symbols, pulse, USF);
% plot(pulse,'.')

r = tx_signal;
r_s = r;
% length(tx_signal) % length(pulse)+#symbols*USF-1 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add impairement simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add random amount of zeros at the begining and end of the signal to simulate begining of transmission at random time
% numZerosStart = randi([1,2000])
% numZerosEnd = randi([1,2000])
% r_s = [zeros(1,numZerosStart),r_s,zeros(1,numZerosEnd)];

%% Add Complex White Gaussian Noise
% complex() makes sure that the signal is considered to be complex so that
% awgn() adds complex noise.
r_s = awgn(complex(tx_signal), SNR, 'measured');

%% Add random initial phase [0,2pi) to simulate time of flight
phi0 = rand*2*pi;
r_s = r_s*exp(1j*phi0);

%% Add constant Doppler shift
fc = 0.1;
dopplerRate = 0.3; %0.5; % 
t = (0:(length(r_s)-1))*Ts;
r_s = r_s.*exp(1j*2*pi*fc*dopplerRate*t);

%% Fixed clock offset simulation
clockOffset = 9; 
USF_s = 10;
r_s = interp(r_s, USF_s);
r_s = downsample(r_s, USF_s, clockOffset);

%% Cut signal at random position
startPos = 1000; %randi([1,length(r_s)/2]);
r_s = r_s(startPos:end);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Due to the Doppler shift, the optimal sampling time of the matched filter has
% changed. We have to estimate the Doppler shift and remove it before
% sampling.
% ...

% Apply matched-filter to obtain sufficient statistics
pulse_MF = conj(fliplr(pulse));
mfOutput = conv(r_s, pulse_MF);

%
% Coarse estimate of the Doppler shift
maxDoppler = 10;
dopplerStep = 0.1;
[doppler_estim, tau_estim] = dopplerCoarseEstimate(mfOutput, fc, Ts, USF, maxDoppler, dopplerStep)

% Remove Doppler shift
t = (0:(length(mfOutput)-1))*Ts;
mfOutput = exp(-1j*2*pi*fc*doppler_estim*t).*mfOutput;

[R,~] = xcorr(mfOutput,upsample(ca,USF));
% R = R(length(mfOutput):end);
% [~,pos] = max(abs(R))
% angle(R(pos))
t3 = (0:length(R)-1);
plot(t3,R)

%%

% Demodulation
% This line  needs to be fixed to account for the amount of actual symbols
symbolEstimates = mfOutput(tau_estim:USF:end-length(pulse)+1);
% symbolEstimates = transpose(symbolEstimates);

dataEstimates = demodulator(symbolEstimates, map);

% Convert 0 to 1 and 1 to 0
% dataEstimates = 

% Plot the constellation

plot(symbolEstimates, 'r*');
hold on;
plot(map, 'kx');
axis([-15,15,-15,15]);
grid on;

size(symbols)
size(symbolEstimates)

% demodulator
% sum(data ~= dataEstimates)
% BER = sum(data ~= dataEstimates)/length(data)
% symbols(1:5)
% symbolEstimates(1:5)
% 
% angle(symbolEstimates(1))

