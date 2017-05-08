function [dataRx,Fs] = simulatorOFDM(dataTx, config)
%SIMULATOROFDM Summary of this function goes here
%   Detailed explanation goes here

Fs = config.Fs;
Ts = 1/Fs; 

channelScaling = 1;

% AWGN
SNR = 30; % dB
awgnEnabled = 1;

% Enable multipath simulation
multipathEnabled = 1;
% amplitudes = randn(1,4); %Gaussian distribution, N(0,1) 
amplitudes = [1.5856, 0.66436, 0.91867, 1.0794];
% amplitudes = [0.47518     0.26142     -1.2723     0.41863]; % with delays
% [0 3 6 12], it makes a frame sync error (1007 instead of 1001), fine
% estimation does NOT solve the error
% amplitudes = [1 0.5];
disp(['Amplitudes of the multipaths: ' num2str(amplitudes)]);
delays = [0 3 6 12]; % in samples (must be less of equal than the cyclic prefix length to avoid ISI)

% Doppler
doppler = 555;

% Clock offset
clockOffsetEnabled = 1;

epsilon = 40; % ppm e.g.2e-6 = 2 ppm

oldVersion = 0;
USFOffsetSimulation = 10;
clockSampleOffset = 5;

% Random cut
randomCutEnabled = 0;
startPos = 123;

dataRx = dataTx.';

if multipathEnabled
    % Multipath channel
    h = createMultipathChannelFilter(amplitudes, delays);
else
    % Simple channel 
    h = 1;
end

% Normalize impulse response
h = h/norm(h);

% Convolve tx_signal with channel impulse response
dataRx = conv(dataRx, h);

% Channel scaling
dataRx = channelScaling*dataRx;

%Add Complex White Gaussian Noise
if awgnEnabled == 1
    dataRx = awgn(complex(dataRx), SNR, 'measured');
end

% Add constant Doppler shift
if doppler ~= 0
    fprintf('Doppler added: %f\n', doppler);
    t = (0:(length(dataRx)-1))*Ts;
    dataRx = dataRx.*exp(1j*2*pi*doppler*t);
else
    disp('No Doppler');
end

% Fixed clock offset simulation
if clockOffsetEnabled == 1
    fprintf('Clock offset added\n');
    if oldVersion
        dataRx = resample(dataRx, USFOffsetSimulation, 1);
        dataRx = downsample(dataRx, USFOffsetSimulation, clockSampleOffset);
    else
        dataRx = resample(dataRx, 1e6+epsilon, 1e6);
    end
else
    disp('No clock offset added');
end

% Cut signal at random position
if randomCutEnabled == 1
    disp('Received signal cut');
    dataRx = dataRx(startPos:end);
else
    disp('Received signal NOT cut');
end

dataRx = dataRx.';

end

