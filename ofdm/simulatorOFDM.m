function [signalRx] = simulatorOFDM(signalTx, Fs, SNR)
%SIMULATOROFDM Allows to simulate different channels and impairments due to
% device imperfections.
% Takes as input:
% - signalTx: the transmitted signal (column vector)
% - Fs: sampling frequency in Hz (scalar)
% Returns:
% - signalRx: the output of the simulated communication channel with
% impairments (column vector)

[signalTx, wasColumnVector] = turnIntoColumnVector(signalTx);

disp('Simulator parameters:');

Ts = 1/Fs;

% Simulate channel scaling
channelScaling = 1;  % {Natural numbers}0
% AWGN
awgnEnabled = 1; % {0,1}
if not(exist('SNR', 'var'))
    SNR = 25; % dB
end

% Multipath channel simulation
multipathEnabled = 1; % {0,1}
amplitudes = randn(1,13); % Gaussian distribution, N(0,1) 
disp(['Amplitudes of the multipaths: ' num2str(amplitudes)]);
delays = [0 1 2 3 4 5 6 7 8 9 10 11 12]; % in samples (must be less of equal than the cyclic prefix length to avoid ISI)

% Frequency offset (aka "Doppler")
frequencyOffset = 456; %Hz

% Sampling Clock Offset (SFO)
samplingClockOffsetEnabled = 0; % {0,1}
epsilon = 40; % ppm e.g. 2 ppm = 2e-6

% Fixed Sample Offset (FSO)
fixedSampleOffsetEnabled = 0; % {0,1}
fixedSampleOffsetUSF = 10;  % {natural numbers}
fixedSampleOffset = 5;  % {0,1,...,fixedSampleOffsetUSF-1}

% Random cut
randomCutEnabled = 0;  % {0,1}
startPos = 123; % {0,1,...,length(signalTx)}

signalRx = signalTx;

% Channel impulse reponse
if multipathEnabled
    % Multipath channel
    h = createMultipathChannelFilter(amplitudes, delays);
else
    % Simple channel 
    h = 1;
end

% Normalize impulse response
h = h/norm(h);

% Convolve signal with channel impulse response
signalRx = conv(signalRx, h);

% Channel scaling
signalRx = channelScaling*signalRx;

%Add Complex White Gaussian Noise
if awgnEnabled
    disp(['SNR: ', num2str(SNR)]);
    signalRx = awgn(complex(signalRx), SNR, 'measured');
end

% Add Frequency Offset (aka Doppler)
if frequencyOffset ~= 0
    fprintf('Frequency Offset added: %f\n', frequencyOffset);
    t = (0:(length(signalRx)-1))*Ts;
    signalRx = signalRx.*exp(1j*2*pi*frequencyOffset*t).';
else
    disp('No Frequency Offset added');
end

% Fixed sample offset simulation
if fixedSampleOffsetEnabled
    fprintf('Fixed Sample Offset added (USF: %d , offset: %d)\n', fixedSampleOffsetUSF, fixedSampleOffset);
    signalRx = resample(signalRx, fixedSampleOffsetUSF, 1);
    signalRx = downsample(signalRx, fixedSampleOffsetUSF, fixedSampleOffset);
else
    disp('No Fixed Sample Offset added');
end
        
% Sampling Frequency Offset
if samplingClockOffsetEnabled
    fprintf('Sampling Clock Offset added (epsilon = %f ppm)\n', epsilon);
    signalRx = resample(signalRx, 1e6+epsilon, 1e6);
else
    disp('No Sampling Frequency Offset added');
end

% Cut signal at random position
if randomCutEnabled
    disp('Received signal cut');
    signalRx = signalRx(startPos:end);
else
    disp('Received signal NOT cut');
end

signalRx = changeVectorBackToRowIfNeeded(signalRx, wasColumnVector);

disp(' ');

end

