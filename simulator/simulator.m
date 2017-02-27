clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulator: ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%...


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controls and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data
load('CAcodes.mat');
ca = satCAcodes(1,:); % Data
N = 40; % Number of repetition of CA code (N > 1)
symbols = repmat(ca,1,N);
numFullCAToKeep = 2; % (> 1) Number of CA to keep (first complete one included)

mode = 1; % 0 means simulation, 1 means USRP

if mode
    % parameters USRP
    fc = 2.4e9; % carrier frequency
    LOO = 100e3; % local oscilator offset
    clockRateTx = 5e6; % main clock
    interpolationTx = 10;
    clockRateRx = clockRateTx; % on the same radio we need to use the same
    decimationRx = interpolationTx; % to get to 1MHz
    clockInputSource = 'Internal';
    OutputDataTypeUSRP = 'double';
    gainTx = 60; %60 loopback 30dB attenuator; % 89 over the air, no attenuator
    gainRx = 40;
    samplesPerFrameRx = 1e5; % max is 375000 in a burst
    loopback = false; % true/false: receive with the same board or not
    burstMode = false; % NOT the same "burst mode" as in the docunentation
    txRxDuration = 3; % in seconds, continous transmission
    
    Ts = interpolationTx/clockRateTx; % 1e-6; % [s] sampling time
    Fs = 1/Ts;
    
    
    noCAtoDropFromStart = 4;
    noCAtoKeep = 1;
    
else
    Fs = 0.5e6;
    Ts = 1/Fs;
end


% Pulse shaping
span = 200;
USF = 5;       % upsampling factor
fprintf('Fs: %f, Ts: %f \n', Fs, Ts);
	 
beta = 0.5;  % roll-off factor
pulse = rcosdesign(beta, span, USF, 'sqrt'); % already normalized to 1
% fvtool(pulse, 'Analysis', 'impulse')   % Visualize the filter

caUpsampled = upsample(ca,USF);

% AWGN
SNR = 30; % Signal-To-Noise ratio in dB
awgnEnabled = 1;

% Doppler
dopplerCorrection = 1;
doppler = 234; % Amount of Doppler for simulation
maxDoppler = 1000; % Absolute value defining the range in which to estimate the Doppler
% dopplerStep = 0.1; % Size of the steps in the range of Doppler
fcorr = 1; % to be removed

% Clock offset
clockOffsetEnabled = 1;
USFOffsetSimulation = 4; % USF used to simulate clock offset
clockSampleOffset = 2; % Offset (in #samples) of the clock
USFClockOffsetCorrection = 10; % USF used to correct for the clock offset

% Random cut
randomCutEnabled = 1;
startPos = length(pulse)+123;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmitter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate symbol-by-symbol pulse train samples
txSignal = symbolsToSamples(symbols, pulse, USF);

rxSignal = txSignal;
rxSignalImpaired = rxSignal;
% length(tx_signal) % length(pulse)+#symbols*USF-1 


if not(mode)
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Impairement simulations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % rxSignalImpaired = [zeros(1,1000), rxSignalImpaired, zeros(1,1000)];    

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
else
    
    %% USRP settings
    
    %SamplesPerFrame = 362;
    %Nframes = floor(length(pufiltered)/SamplesPerFrame)
    
    pufiltered = rxSignalImpaired; 
    Ts = interpolationTx/clockRateTx;
    
    % find connected radios
    connectedRadios = findsdru
    
    usrpBoardPlatform = 'B200'; % connectedRadios.Platform;
    usrpBoardSerialNumTx = '30C51BC'; % connectedRadios.SerialNum;
    if loopback
        usrpBoardSerialNumRx = '30C51BC'; %'30C5426'; % usrpBoardSerialNumTx;
    else
        usrpBoardSerialNumRx = '30C5426'; %'30C5426'; % usrpBoardSerialNumTx;
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
        'UnderrunOutputPort', true)
    
    if burstMode
        tx.EnableBurstMode = true;
        no_frames = 1; % floor(length(pufiltered)/samplesPerFrameRx) % for other TxRx option, does not seem to work
        tx.NumFramesInBurst = no_frames;
    end
    
    TxInfo = info(tx)
    
    % configure the Rx object
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
        'OutputDataType', OutputDataTypeUSRP,...
        'OverrunOutputPort', true)
    
    if burstMode
        rx.EnableBurstMode = true;
        rx.NumFramesInBurst = no_frames;
    end
    
    RxInfo = info(rx)
    
    fprintf(1,'Length of the Tx data: %d. \n', length(pufiltered));
    
    % here is with bursts
    % dataRx = zeros(1, length(pufiltered));
    % loop over the number of frames
    % for indexFrame = 1:Nframes
    %     % transmit
    %     dataToTx = pufiltered((indexFrame-1)*SamplesPerFrame+1:indexFrame*SamplesPerFrame).';
    %     %dataTxUSRP = tx(dataToTx);
    %     step(tx, dataToTx);
    %     % start receiver
    %     %dataRxUSRP = rx();
    %     dataRxUSRP = step(rx);
    %     dataRx((indexFrame-1)*SamplesPerFrame+1:indexFrame*SamplesPerFrame) = dataRxUSRP;
    % end
    
    figure;
    plot(pufiltered); grid on; title('Tx data, after Pulse Shaping');
    
    % try two options for Tx/Rx
    if burstMode
        %% old way, not very good
        % here is continous: we loop to have enough data at Rx
        for k = 1:6
            dataTxUSRP  = tx(pufiltered.');
            [dataRxUSRP,len,lostSamples] = rx();
            fprintf(1, 'Receive frame %i in the burst, lost samples at Rx: %i. \n', k, lostSamples);
        end
        dataRx = dataRxUSRP.';
        
    else
        %% new way from Ethem Sozer @mathworks
        t = 0;
        tFrame = rx.SamplesPerFrame/(rx.MasterClockRate/rx.DecimationFactor);
        while t < txRxDuration
            [dataRxUSRP, len] = rx();
            if len > 0
                dataRx = dataRxUSRP.';
                dataTxUSRP = tx(pufiltered.');
                t = t + tFrame;
            end
        end
        
    end
    
    %% yet another way as in usrpReceiveSamples.m
    %     totalSamples = no_frames*samplesPerFrameRx;
    %
    %     for k = 1:3
    %         toto = zeros(totalSamples, 1);
    %         samplesLost = 0;
    %
    %         for l = 1:no_frames
    %             dataTxUSRP  = tx(pufiltered(1+(l-1)*samplesPerFrameRx:l*samplesPerFrameRx).');
    %             [dataRxUSRP,len,lostSamples] = rx();
    %             samplesLost = samplesLost + lostSamples;
    %             %fprintf(1, 'File %i, frame %i, lost samples at Rx: %i. \n', k, l, lostSamples);
    %             toto(1+(l-1)*samplesPerFrameRx:l*samplesPerFrameRx) = dataRxUSRP;
    %         end
    %     end
    %
    %     fprintf(1, 'Lost samples at Rx: %i. \n', samplesLost);
    %
    %     fileName = 'dataRxTmp.mat';
    %     save(fileName, '-v7.3', 'toto');
    %     dataRx = toto.';
    %
    
    
    %%
    % release the objects
    release(tx)
    release(rx)
    
    fprintf(1, 'Length of Tx data for USRP: %d. \n', length(dataTxUSRP));
    fprintf(1, 'Length of Rx data for USRP: %d. \n', length(dataRxUSRP));
    
    % dataRxUSRP = dataRxUSRP - mean(dataRxUSRP); % remove DC offset
    % dataRxUSRP = dataRxUSRP/var(dataRxUSRP); % set to var 1
    
    scatterplot(dataRx);
    grid on; title('Received Sequence');
    
    % drop some data
    lengthCA_SPS = length(caUpsampled);
    startData = noCAtoDropFromStart*lengthCA_SPS;
    stopData = startData + (noCAtoKeep+2)*lengthCA_SPS;
    dataRx1 = dataRx(startData:stopData);
    
    figure;
    hax=axes;
    plot(abs(dataRx)); grid on; title('abs(Received Sequence): vertical bars show selected data')
    hold on;
    VL1 = startData;
    VL2 = stopData;
    line([VL1 VL1],get(hax,'YLim'), 'Color', [1 0 0]);
    line([VL2 VL2],get(hax,'YLim'), 'Color', [1 0 0]);
    
    %dataRx = dataRx - mean(dataRx);
    
    rxSignalImpaired = dataRx1;
    
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matched-Filter
matchedFilter = conj(fliplr(pulse));
MFOutput = conv(rxSignalImpaired, matchedFilter);

scatterplot(MFOutput);
title('Matched-filter output');

figure;
plot(abs(MFOutput))
title('Magnitude of matched-filter output');
xlabel('Sample number');
ylabel('Magnitude');

% Upsample Rx signal to correct for clock offset
caUpsampledClockOffset = upsample(caUpsampled, USFClockOffsetCorrection);
MFOutputUpsampled = resample(MFOutput, USFClockOffsetCorrection, 1);

% Doppler shift estimation and correction
fprintf('\n Doppler shift estimation and correction');
signalCorrected = dopplerEstimationAndCorrection(MFOutput, caUpsampled, Ts, maxDoppler, fcorr, dopplerCorrection);

% Doppler shift estimation and correction from upsampled signal (to account
% for the clock offset as well)
fprintf('\n Doppler shift estimation and correction for upsampled signal');
signalUpsampledCorrected = dopplerEstimationAndCorrection(MFOutputUpsampled, caUpsampledClockOffset, Ts/USFClockOffsetCorrection, maxDoppler, fcorr, dopplerCorrection);

% Sample the matched filter to get the symbols estimates
symbolsEstimates = signalCorrected(1:USF:length(caUpsampled)*numFullCAToKeep-1);
symbolsUpsampledEstimates = signalUpsampledCorrected(1:USF*USFClockOffsetCorrection:length(caUpsampledClockOffset)*numFullCAToKeep-1);

% Plot the constellation
for k = 1:numFullCAToKeep
    scatterplot(symbolsEstimates(1:length(ca)*k));
    title(['BPSK constellation: first ' num2str(k) ' CA code(s)']);
    hold on;
    plot([-1 1],[0 0], 'rx');
    grid on;
    
    scatterplot(symbolsUpsampledEstimates(1:length(ca)*k));
    title(['(UpS) BPSK constellation: first ' num2str(k) ' CA code(s)']);
    hold on;
    plot([-1 1],[0 0], 'rx');
    grid on;
end

map = pskMap(2);
% size(symbols(1:length(ca)*numFullCAToKeep))
% size(symbolsEstimates)
symbolsDec = demodulator(symbols(1:length(ca)*numFullCAToKeep),map);
symbolsEstimatesDec = demodulator(symbolsEstimates,map);
BER = sum(symbolsDec ~= symbolsEstimatesDec)/length(symbolsDec);

fprintf('\n BER: %f \n', BER);

