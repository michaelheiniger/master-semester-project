clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit and receive bits using OFDM modulation
% (See the bottom of this file for acronyms and abbreviations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Instance started on %s \n\n', datestr(now))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controls and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Possible modes: 
% - 'simulation'        can simulate different impairments: AWGN, CFO, SFO,

mode = 'simulation';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System configuration 
sc.Fs = 0.5e6; % [Hz], sampling frequency
sc.M = 4; % Size of M-QAM constellation
sc.map = qammap(sc.M); % M-QAM constellation vector

sc.numTotalCarriers = 64; % Total number of subcarriers, guard bands included (FFT size)
sc.NFFT = sc.numTotalCarriers; % FFT size
sc.numPilots = 2; % Number of subcarriers used as pilots (NOTE: system won't adapt to change)
sc.CPLength = 16; % Length of the cyclic prefix of the pilot and regular OFDM symbols (16 in IEEE 802.11a)
sc.twoLtsCpLength = 32; % Length of the cyclic prefix of the two LTSs of the preamble (32 in IEEE 802.11a)
sc.numOFDMSymbolsPerFrame = 200; % Excluding preamble, must be >= 2

marginSubCarriers = 0.17*sc.numTotalCarriers; % Ratio of subcarriers that are considered to be at the margin of the spectrum used
sc.numZerosTop = ceil(marginSubCarriers/2); % Amount of top (negative freq) subcarriers used as guard bands
sc.numZerosBottom = floor(marginSubCarriers/2); % Amount of bottom (positive freq) subcarriers used as guard bands
% config.numZerosTop = 1; > 0 since channel estimation is done on guard bands
% config.numZerosBottom = 1; > 0 since channel estimation is done on guard bands
sc.zeroFreqSubcarrier = 1; % 0 if used for data, 1 means that it is set to zero

sc.numUsedCarriers = sc.numTotalCarriers - sc.numZerosTop - sc.numZerosBottom - sc.zeroFreqSubcarrier; % number of subcarriers used to carry symbols (information or pilot)
sc.numDataCarriers = sc.numUsedCarriers - sc.numPilots; % number of subcarriers used to carry information symbols (i.e. not pilots)
sc.numBitsForPayloadSize = 16; % first bits of signal field, must be < numDataCarriers

[stsTime, ~] = getSTS();
[ltsTime, ~] = getLTS();
sc.ofdmFrameLength = 10*length(stsTime) + sc.twoLtsCpLength + 2*length(ltsTime) + sc.numOFDMSymbolsPerFrame*(sc.CPLength+sc.NFFT);
sc.usrpFrameLength = 5e4;
sc % Display config

% The noise variance is estimated using the guard bands so there must be at
% least one.
if sc.numZerosBottom+sc.numZerosTop == 0
   error('There is no zero subcarriers: cannot measure the noise variance using MMSE !'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Receivers configuration
% Several receivers can be defined for comparison
% Methods for timing and frequency offset estimation are:
% - stsLtsOfdmDemod: use IEEE 802.11a preamble and demodulates preamble
%   symbols to compute estimates
% - ideal use manually chosen timing and CFO estimates.

% Configuration of OFDM receiver 1
rc1.timingAndFrequencyOffsetMethod = 'stsLtsOfdmDemod';
rc1.timingThreshold = 0.3
rc1.cfoCorrection = 1; % 1 if CFO should be corrected
rc1.cfoTracking = 1; % 1 if residual CFO should be tracked
rc1.sfoCorrection = 0; % 1 if SFO should be corrected
rc1.equalization = 1; % 1 if channel equalization should be performed
rc1.timingOffset = 0; % add offset to timing estimate, in number of samples.

% Configuration of OFDM receiver 2
rc2.timingAndFrequencyOffsetMethod = 'ideal';
rc2.manualCFO = 456; % Set the timing sample regardless of the actual estimate
rc2.cfoCorrection = 1; % 1 if CFO should be corrected
rc2.cfoTracking = 1; % 1 if residual CFO should be tracked
rc2.sfoCorrection = 0; % 1 if SFO should be corrected
rc2.equalization = 1; % 1 if channel equalization should be performed
rc2.timingOffset = 0; % add offset to timing estimate, in number of samples.

numRuns = 1; % Different random data from one run to the other
numRunsPerSnr = 25; % Different noise and channel from one run to the other
snrValues = 20:1:55; % Different SNR values
% numRuns = 1; % Different random data from one run to the other
% numRunsPerSnr = 1; % Different noise and channel from one run to the other
% snrValues = 20; % Different SNR values
receiver1Results = zeros(numRuns, length(snrValues), numRunsPerSnr);
receiver1Timings = zeros(numRuns, length(snrValues), numRunsPerSnr);
receiver2Results = zeros(numRuns, length(snrValues), numRunsPerSnr);

for run = 1:numRuns
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Transmitter
    % Fetch the bits to send and builds the OFDM frames
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear signalTx;
    
    if isInstanceTransmitter(mode)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % OFDM transmitter
        % Build the OFDM frames and serialize them into signalTx
        % dataFrame is the frame containing the subcarriers used to carry
        % symbols (information or pilot) and is returned only for comparison
        % with received version (see OFDMReceiver)
        bitsToSend = [];
        [signalTx, dataFrame, bitsToSendWithPadding, infoSymbols] = OFDMTransmitter(bitsToSend, sc);

    else
        % If the current MATLAB instance if NOT a transmitter, the signal to
        % transmit is empty (necessary for USRP code)
        signalTx = [];
    end
    
    for snrIndex = 1:length(snrValues)
        for snrRun = 1:numRunsPerSnr
            clear signalRx;
            clear infoSymbolsRx1;
            clear bitsRx1;
            clear infoSymbolsRx2;
            clear bitsRx2;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Transmission / Reception of the signal over the channel
            % Simulator or USRP
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            plotSignalMagnitude(signalTx, 'Samples', 'Absolute value of transmitted signal')

            % Simulate the impairements of the channel on the signal
            [signalTxImpaired, delayLargestAmplitude] = simulatorOFDM(signalTx, sc.Fs);
            
            % Set the timing manually for the ideal receiver: the coarse
            % frame sync will fix on the strongest copy of the signal. The
            % delay of this copy is known to the ideal receiver so that it
            % can adapt and always use the true OFDM frame beginning
            rc2.manualTiming = 1 + sc.CPLength - delayLargestAmplitude;

            % Place the OFDM frames at the beginning of the received signal if
            % deterministic = 1, else the frames are randomly placed
            deterministic = 1;
            [coarseFrameRx, ~, ~] = usrpSimulator(signalTxImpaired, sc.usrpFrameLength, sc.ofdmFrameLength, sc.CPLength, deterministic);    
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % OFDM Receivers
            % Feed the different OFDM receivers with the received signal
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isInstanceReceiver(mode)
                
                % Receiver 1
                [bitsRx1, infoSymbolsRx1, ~, timingEst1] = OFDMReceiver(sc, rc1, coarseFrameRx, dataFrame);
                
                 % Compute Bit and Symbol Error Rates of all receivers
                [BER1, SER1] = symbolsAndBitsStats(infoSymbols, infoSymbolsRx1, bitsToSendWithPadding, bitsRx1, sc.map)

                % Save results for receiver 1
                receiver1Results(run, snrIndex, snrRun) = SER1;
                receiver1Timings(run, snrIndex, snrRun) = timingEst1;
                
                % Receiver 2
                [bitsRx2, infoSymbolsRx2, ~, ~] = OFDMReceiver(sc, rc2, coarseFrameRx, dataFrame);

                [BER2, SER2] = symbolsAndBitsStats(infoSymbols, infoSymbolsRx2, bitsToSendWithPadding, bitsRx2, sc.map)

                % Save results for receiver 2
                receiver2Results(run, snrIndex, snrRun) = SER2;
                    
            end
            if numRunsPerSnr > 1
                close all;
            end
        end
    end
    
end

%% Save results
resultFileReceiver1 = 'receiver1Results_batch3.mat';
timingResultFileReceiver1 = 'timingResultFileReceiver1_batch3.mat';
resultFileReceiver2 = 'receiver2Results_batch3.mat';

save(resultFileReceiver1, 'receiver1Results');
save(resultFileReceiver2, 'receiver2Results');
save(timingResultFileReceiver1, 'receiver1Timings');


%% Compute mean over main runs
meanOverData1 = mean(receiver1Results,1);
meanOverData2 = mean(receiver2Results,1);

% Compute variance over SNR runs (i.e. variance of noise and multipath)
varResultsReceiver1 = var(meanOverData1,0,3);
varResultsReceiver2 = var(meanOverData2,0,3);

% Compute mean over SNR runs
meanOverDataThenSnr1 = mean(meanOverData1, 3);
meanOverDataThenSnr2 = mean(meanOverData2, 3);

%%
figure;
errorbar(snrValues,meanOverDataThenSnr1,varResultsReceiver1)
hold on;
errorbar(snrValues,meanOverDataThenSnr2,varResultsReceiver2)
% axis([15,55]);
xlabel('SER');
ylabel('SNR [dB]');
legend('Actual Rx','Ideal Rx');
title('SER vs SNR');

fprintf('Instance terminated on %s \n\n',datestr(now));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acronyms and abbreviations:
% AWGN: Additive White Gaussian Noise
% CFO: Carrier Frequency Offset
% SFO: Sampling Frequency Offset
% USRP: Universal Software Radio Peripheral
% Tx: Transmitter
% Rx: Receiver
% CP: Cyclic Prefix
% STS: Short Training Sequence
% LTS: Long Training Sequence
% CA: Coarse Acquisition codes (in reference to GPS C/A codes)
% USF: Up-Sampling Factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%