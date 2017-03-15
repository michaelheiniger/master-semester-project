clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to send and receive data through simulator or USRP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Instance started on %s \n\n',datestr(now))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controls and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Possible modes: 
% - 'simulation':       simulates impairements such as AWGN, Doppler, clock offset,... 
% - 'loopback':         one USRP board receives its own transmission 
%                       --> Tx and Rx ports need to be linked by cable 
% - 'oneBoardTx':       one USRP board transmits 
% - 'oneBoardRx':       one USRP board receives
% - 'twoBoardsRxTx':    two USRP boards on the same computer: one transmits, the
%                       other receives
% USRP boards are always used in half-duplex mode

mode = 'simulation'; 
fprintf('Mode of operation: %s \n\n', mode)

burstMode = 0 ;
fprintf('Burst mode enabled: %s \n\n', num2str(burstMode))

% Source data
load('CAcodes.mat');
ca = satCAcodes(1,:);
N = 1000; % Number of repetition of CA code (N > 1)
symbols = repmat(ca,1,N);

noCAtoDropFromStart = 40;
noCAtoKeep = 50; % (> 1) Number of CA to keep (first complete one included)
idCaToShow = [1 5 10 20 30 40];

samplesPerFrame = 1e5; % max is 375000 in a burst

% Pulse shaping
span = 200; 
USF = 5; % upsampling factor
	 
beta = 0.8;  % roll-off factor
pulse = rcosdesign(beta, span, USF, 'sqrt'); % already normalized to 1
% fvtool(pulse, 'Analysis', 'impulse')   % Visualize the filter

caUpsampled = upsample(ca,USF);

% Doppler
dopplerCorrection = 1;
maxDoppler = 2000; % Absolute value defining the range in which to estimate the Doppler
    
% Clock offset
USFClockOffsetCorrection = 10; % USF used to correct for the clock offset

% Generate symbol-by-symbol pulse train samples
dataTx = symbolsToSamples(symbols, pulse, USF);

% Pad with zeros to make an integer number of frames
dataTx = padLastFrameWithZeros(dataTx, samplesPerFrame);

noFrames = length(dataTx)/samplesPerFrame;

if strcmp(mode, 'loopback') || strcmp(mode, 'oneBoardTx') ||  strcmp(mode, 'twoBoardsRxTx')
    fprintf('Transmit: %s frames for %s CA codes \n', num2str(noFrames), num2str(N));
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmitter / Receiver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;
% plot(dataTx);
% title('Signal to transmit')

switch mode
    case {'oneBoardRx','oneBoardTx','twoBoardsRxTx','loopback'}
        [dataRx, Fs] = rxTxUSRP(dataTx, samplesPerFrame, mode);
    case 'simulation'
        [dataRx, Fs] = simulatorPS(dataTx, length(pulse));
    otherwise
        error('Unknown mode of operation.')
end

Ts = 1/Fs;

if strcmp(mode, 'simulation') || strcmp(mode, 'loopback') || strcmp(mode, 'oneBoardRx') ||  strcmp(mode, 'twoBoardsRxTx')
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Signal equalization and decoding
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf(1, 'Length of Rx data: %d. \n', length(dataRx));

    % scatterplot(dataRx);
    % grid on; 
    % title('Received Sequence');

    % drop some data
    lengthCA_SPS = length(caUpsampled);
    startData = noCAtoDropFromStart*lengthCA_SPS;
    stopData = length(dataRx); %startData + noCAtoKeep*lengthCA_SPS;
    dataRx1 = dataRx(startData:stopData);

    rxSignalImpaired = dataRx1;

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
    MFOutputUpsampled = resample(MFOutput, USFClockOffsetCorrection, 1);


    % Doppler shift estimation and correction
    fprintf('\nDoppler shift estimation and correction');
    [signalCorrected, tauEst] = dopplerEstimationAndCorrection(MFOutput, caUpsampled, Ts, maxDoppler, dopplerCorrection);

    % Doppler shift estimation and correction from upsampled signal (to account
    % for the clock offset as well)
    fprintf('\n Doppler shift estimation and correction for upsampled signal \n');
    [signalUpsampledCorrected, ~] = dopplerEstimationAndCorrection(MFOutputUpsampled, caUpsampledClockOffset, Ts/USFClockOffsetCorrection, maxDoppler, dopplerCorrection);

    figure;
    hax=axes;
    plot(abs(dataRx)); grid on; title('abs(Received Sequence): red bars show selected data and green bar shows Tau')
    hold on;
    VL1 = startData;
    VL2 = stopData;
    VL3 = startData + tauEst - 1;
    VL4 = startData + tauEst + length(caUpsampled)*noCAtoKeep - 1;
    line([VL1 VL1],get(hax,'YLim'), 'Color', [1 0 0]);
    line([VL2 VL2],get(hax,'YLim'), 'Color', [1 0 0]);
    line([VL3 VL3],get(hax,'YLim'), 'Color', [0 1 0]);
    line([VL4 VL4],get(hax,'YLim'), 'Color', [0 1 0]);


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

    fprintf('\nStatistics \n');

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
end

fprintf('Instance terminated on %s \n\n',datestr(now))


