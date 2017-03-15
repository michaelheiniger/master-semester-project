clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulator: Script to send and receive data through simulator or USRP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Instance started on %s \n\n',datestr(now))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controls and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Switch between simulator and USRP
mode = 1; % 0 for simulator, 1 for USRP

% Source data
load('CAcodes.mat');
ca = satCAcodes(1,:);
N = 1000; % Number of repetition of CA code (N > 1)
symbols = repmat(ca,1,N);

noCAtoDropFromStart = 40;
noCAtoKeep = 50; % (> 1) Number of CA to keep (first complete one included)
idCaToShow = [1 5 10 20 30 40];

% Pulse shaping
span = 200;
USF = 5;       % upsampling factor
	 
beta = 0.8;  % roll-off factor
pulse = rcosdesign(beta, span, USF, 'sqrt'); % already normalized to 1
% fvtool(pulse, 'Analysis', 'impulse')   % Visualize the filter

caUpsampled = upsample(ca,USF);

% Doppler
dopplerCorrection = 1;
maxDoppler = 2000; % Absolute value defining the range in which to estimate the Doppler
    
% Clock offset
USFClockOffsetCorrection = 10; % USF used to correct for the clock offset

if mode
    % USRP parameters
    fc = 2.4e9; % carrier frequency
    LOO = 100e3; % local oscillator offset
    clockRateTx = 5e6; % main clock (Sampling rate of the digital signal sent to ADC) (5e6 to 56e6)
    interpolationTx = 10; % Interpolation factor from host signal to USRP signal (e.g Fs = 0.5 Mhz => ClockRate = 5 Mhz)
    clockRateRx = clockRateTx; % on the same radio we need to use the same
    decimationRx = interpolationTx; % to get to 1MHz
    clockInputSource = 'Internal';
    outputDataTypeUSRP = 'double'; % Difference between transport and output data type ?
    gainTx = 60; %60 loopback 30dB attenuator; % 89 over the air, no attenuator
    gainRx = 40;
    samplesPerFrameRx = 1e5;%205600; %205600 = signal of 40 CA codes %1e5; % max is 375000 in a burst
    loopback = true; % true/false: receive with the same board or not
    burstMode = false;
    sendMode = false;
    
    Ts = interpolationTx/clockRateTx; % [s] sampling time 
    Fs = 1/Ts;
    
%     txRxDuration = 3; % in seconds, continous transmission
    
else
    % Simulator parameters
    Fs = 0.5e6;
    Ts = 1/Fs;
    
    % AWGN
    SNR = 30; % Signal-To-Noise ratio in dB
    awgnEnabled = 1;
    
    % Doppler
    doppler = 234; % Amount of Doppler for simulation
    
    % Clock offset
    clockOffsetEnabled = 1;
    USFOffsetSimulation = 4; % USF used to simulate clock offset
    clockSampleOffset = 2; % Offset (in #samples) of the clock
    
    % Random cut
    randomCutEnabled = 1;
    startPos = length(pulse)+123;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmitter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate symbol-by-symbol pulse train samples
txSignal = symbolsToSamples(symbols, pulse, USF);

% figure;
% plot(txSignal);
% title('Signal to transmit')
if mode % mode==1: USRP, mode==0: Simulation
    %% USRP settings
    
    pufiltered = txSignal;
    
    % find connected radios
    connectedRadios = findsdru
    
    usrpBoardPlatform = 'B200'; % connectedRadios.Platform;
    usrpBoardSerialNumTx = '30C5426'; % connectedRadios.SerialNum;
    if loopback
        usrpBoardSerialNumRx = '30C5426'; % '30C51BC','30C5426'; % usrpBoardSerialNumTx;
    else
        usrpBoardSerialNumRx = '30C51BC'; %'30C5426'; % usrpBoardSerialNumTx;
    end
    
        
    % configure the Tx object
    tx = comm.SDRuTransmitter(...
        'Platform',usrpBoardPlatform, ...
        'CenterFrequency',fc, ...
        'LocalOscillatorOffset', LOO,...
        'InterpolationFactor',interpolationTx, ...
        'MasterClockRate', clockRateTx,...
        'SerialNum', usrpBoardSerialNumTx,...
        'Gain', gainTx,...
        'ClockSource', clockInputSource,...
        'UnderrunOutputPort', true); % Signal when underrun occurs
    
    if burstMode
        tx.EnableBurstMode = true;
        noFrames = 20; % floor(length(pufiltered)/samplesPerFrameRx) % for other TxRx option, does not seem to work
        tx.NumFramesInBurst = noFrames;
    end
    
    TxInfo = info(tx);
    
    % Configure the Rx object
    rx = comm.SDRuReceiver(...
        'Platform',usrpBoardPlatform, ...
        'CenterFrequency',fc, ...
        'LocalOscillatorOffset', LOO,...
        'MasterClockRate', clockRateRx,...
        'SerialNum', usrpBoardSerialNumRx,...
        'DecimationFactor', decimationRx,...
        'Gain', gainRx,...
        'SamplesPerFrame', samplesPerFrameRx,...
        'ClockSource', clockInputSource,...
        'OutputDataType', outputDataTypeUSRP,...
        'OverrunOutputPort', true); % Signal when overrun occurs
    
    if burstMode
        rx.EnableBurstMode = true;
        rx.NumFramesInBurst = noFrames;
    end
    
    RxInfo = info(rx);
    
    fprintf(1,'Length of the Tx data: %d. \n', length(pufiltered));
    
    
 

        noFrames = floor(length(pufiltered)/rx.SamplesPerFrame); %floor(txRxDuration/tFrame);
        dataRx = zeros(1, noFrames*rx.SamplesPerFrame);
        l = 1; w = 1;
        while l <= noFrames
            [dataRxUSRP, len, overrun] = rx();

            if len > 0
                dataRx(1+(l-1)*rx.SamplesPerFrame:l*rx.SamplesPerFrame) = dataRxUSRP.';
                dataTxUSRP = tx(pufiltered(1+(l-1)*rx.SamplesPerFrame:l*rx.SamplesPerFrame).');
                l = l+1;
            end
        end
        
    % release the objects
    release(tx)
    release(rx)
    
    fprintf(1, 'Length of Rx data for USRP: %d. \n', length(dataRx));
    

%     scatterplot(dataRx);
%     grid on; 
%     title('Received Sequence');
    
    % drop some data
    lengthCA_SPS = length(caUpsampled);
    startData = noCAtoDropFromStart*lengthCA_SPS;
    stopData = length(dataRx); %startData + noCAtoKeep*lengthCA_SPS;
    dataRx1 = dataRx(startData:stopData);
    
    figure;
    hax=axes;
    plot(abs(dataRx)); grid on; title('abs(Received Sequence): vertical bars show selected data')
    hold on;
    VL1 = startData;
    VL2 = stopData;
    line([VL1 VL1],get(hax,'YLim'), 'Color', [1 0 0]);
    line([VL2 VL2],get(hax,'YLim'), 'Color', [1 0 0]);

    rxSignalImpaired = dataRx1;

else % SIMULATION
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Impairement simulations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dataRx1 = txSignal;
    rxSignalImpaired = txSignal;

    figure;
    plot(abs(dataRx1));
    
    %Add Complex White Gaussian Noise
    if awgnEnabled == 1
        rxSignalImpaired = awgn(complex(rxSignalImpaired), SNR, 'measured');
    end

    % Add constant Doppler shift
    if doppler ~= 0
        fprintf('Doppler added: %f\n', doppler);
        t = (0:(length(rxSignalImpaired)-1))*Ts;
        rxSignalImpaired = rxSignalImpaired.*exp(1j*2*pi*doppler*t);
    else
       disp('No Doppler'); 
    end

    % Fixed clock offset simulation
    if clockOffsetEnabled == 1
        disp('Clock offset added'); 
        rxSignalImpaired = resample(rxSignalImpaired, USFOffsetSimulation, 1);
        rxSignalImpaired = downsample(rxSignalImpaired, USFOffsetSimulation, clockSampleOffset);
    else
        disp('No clock offset added'); 
    end

    % Cut signal at random position
    if randomCutEnabled == 1
        disp('Received signal cut'); 
        rxSignalImpaired = rxSignalImpaired(startPos:end);
    else
        disp('Received signal NOT cut'); 
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Apply Matched-Filter
matchedFilter = conj(fliplr(pulse));
MFOutput = conv(rxSignalImpaired, matchedFilter);

scatterplot(MFOutput);
title('Matched-filter output');

% figure;
% plot(abs(MFOutput))
% title('Magnitude of matched-filter output');
% xlabel('Sample number');
% ylabel('Magnitude');

% Upsample Rx signal to correct for clock offset
caUpsampledClockOffset = upsample(caUpsampled, USFClockOffsetCorrection);
% caUpsampledClockOffset = hAGC(upsample(caUpsampled.', USFClockOffsetCorrection)).';
MFOutputUpsampled = resample(MFOutput, USFClockOffsetCorrection, 1);

% Doppler shift estimation and correction
fprintf('\n Doppler shift estimation and correction');

[signalCorrected, tauEst] = dopplerEstimationAndCorrection(MFOutput, caUpsampled, Ts, maxDoppler, dopplerCorrection);
% [signalCorrected, tauEst] = dopplerEstimationAndCorrection(MFOutput, hAGC(caUpsampled.').', Ts, maxDoppler, dopplerCorrection);

figure;
hax=axes;
plot(abs(dataRx1)); grid on; title('abs(Received Sequence): vertical bars show selected data')
hold on;
VL1 = tauEst;
line([VL1 VL1],get(hax,'YLim'), 'Color', [0 1 0]);

% Doppler shift estimation and correction from upsampled signal (to account
% for the clock offset as well)
fprintf('\n Doppler shift estimation and correction for upsampled signal \n');
[signalUpsampledCorrected, ~] = dopplerEstimationAndCorrection(MFOutputUpsampled, caUpsampledClockOffset, Ts/USFClockOffsetCorrection, maxDoppler, dopplerCorrection);

% Sample the matched filter to get the symbols estimates
symbolsEstimates = signalCorrected(1:USF:length(caUpsampled)*noCAtoKeep-1);
symbolsUpsampledEstimates = signalUpsampledCorrected(1:USF*USFClockOffsetCorrection:length(caUpsampledClockOffset)*noCAtoKeep-1);

% Plot the constellation
for k = idCaToShow
%     scatterplot(symbolsEstimates(1:length(ca)*k));
%     title(['BPSK constellation: first ' num2str(k) ' CA code(s)']);
%     hold on;
%     plot([-1 1],[0 0], 'rx');
%     grid on;
    
    scatterplot(symbolsUpsampledEstimates(1:length(ca)*k));
    title(['(UpS) BPSK constellation: first ' num2str(k) ' CA code(s)']);
    hold on;
    plot([-1 1],[0 0], 'rx');
    grid on;
end

fprintf('\n Statistics \n');

% TO CORRECT !!!!
% length(pulse)+length(symbols)*USF-1 
% bitRate = Fs/length(pulse); % bits/s
% fprintf('Bit rate: %f bits/s \n', bitRate);
bandwidth = Fs*(1+beta); % Hz
fprintf('Bandwidth used: %f Mhz\n', bandwidth*1e-6);

map = pskMap(2);
% size(symbols(1:length(ca)*numCAtoKeep))
% size(symbolsEstimates)
symbolsDec = demodulator(symbols(1:length(ca)*noCAtoKeep),map);
symbolsEstimatesDec = demodulator(symbolsEstimates,map);
symbolsUpsampledEstimatesDec = demodulator(symbolsUpsampledEstimates,map);

fprintf('Number of symbols decoded: %d \n', length(ca)*noCAtoKeep);
BER = sum(symbolsDec ~= symbolsEstimatesDec)/length(symbolsDec);
BERUpsampled = sum(symbolsDec ~= symbolsUpsampledEstimatesDec)/length(symbolsDec);

fprintf('BER: %f \n', BER);
fprintf('BER (Upsampled): %f \n', BERUpsampled);

