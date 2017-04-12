function [dataRx,Fs] = simulatorOFDM(dataTx)
%SIMULATOROFDM Summary of this function goes here
%   Detailed explanation goes here

global decoder;

Fs = 0.5e6;
Ts = 1/Fs; 

% AWGN
SNR = 30;
awgnEnabled = 1;

% Doppler
doppler = 345

% Clock offset
clockOffsetEnabled = 1;

epsilon = 100; % ppm e.g.2e-6 = 2 ppm

oldVersion = 0;
USFOffsetSimulation = 10;
clockSampleOffset = 5;

% Random cut
randomCutEnabled = 0;
startPos = 123;

dataRx = dataTx.';

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
    fprintf('Clock offset added (%d for USF %d)\n', clockSampleOffset, USFOffsetSimulation);
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

