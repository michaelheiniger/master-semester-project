function [dataRx,Fs] = simulatorPulseShaping(dataTx, lengthPulse)
%SIMULATORPULSESHAPING Summary of this function goes here
%   Detailed explanation goes here

Fs = 0.5e6;
Ts = 1/Fs;

% AWGN
SNR = 30;
awgnEnabled = 1;

% Doppler
doppler = 234;

% Clock offset
clockOffsetEnabled = 1;
USFOffsetSimulation = 4;
clockSampleOffset = 2;

% Random cut
randomCutEnabled = 1;
startPos = lengthPulse+123;

dataRx = dataTx;

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
    disp('Clock offset added');
    dataRx = resample(dataRx, USFOffsetSimulation, 1);
    dataRx = downsample(dataRx, USFOffsetSimulation, clockSampleOffset);
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


end

