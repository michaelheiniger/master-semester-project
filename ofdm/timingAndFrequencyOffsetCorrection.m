function [frameCorrected] = timingAndFrequencyOffsetCorrection(systemConfig, receiverConfig, signal, stsTime, ltsTime, ca)
%TIMINGANDFREQUENCYOFFSET Summary of this function goes here
%   Detailed explanation goes here

sc = systemConfig;
rc = receiverConfig;

Ts = 1/sc.Fs;

% Timing and Frequency Offset Estimation
% - timingEst is the estimation of the index of the first sample of the frame
% - frequencyOffsetEst is the frequency offset the time-dependent complex exponential introducing
%   phase-shift drifting over time due to carrier frequency offset between
%   transmitter and receiver.
% - intialPhaseEst is the constant phase-shift part of the frequency offset
%   The point of reference for the inital phase is the start of the frame
% Here are several methods to compare, see the MATLAB functions header for
% more details
switch rc.timingAndFrequencyOffsetMethod
    case 'stsLtsOfdmDemod'
        [timingEst, frequencyOffsetEst, initialPhaseEst] = timingAndFrequencyOffsetEstimation1(sc, signal, stsTime, ltsTime, ca);
    case 'stsLtsTimeDomain'
        [timingEst, frequencyOffsetEst, initialPhaseEst] = timingAndFrequencyOffsetEstimation3(sc, signal, stsTime, ltsTime, ca);
    case 'caTimeDomain'
        [timingEst, frequencyOffsetEst, initialPhaseEst] = timingAndFrequencyOffsetEstimation2(sc, rc, transpose(signal), ca);
    case 'ideal'
        timingEst = rc.manualTiming;
        frequencyOffsetEst = rc.manualCFO;
        initialPhaseEst = angle(exp(1j*2*pi*frequencyOffsetEst*(timingEst-1)*Ts));
    otherwise
        error('Timing and frequency offset correction method unknown.');
end

disp(['Timing: ', num2str(timingEst)]);
disp(['Total CFO: ', num2str(frequencyOffsetEst)]);
disp(['Initial phase: ', num2str(initialPhaseEst)]);

% Allows to introduce manually a timing offset (see Receiver configuration)
timingEst = timingEst + rc.timingOffset;

lengthStsTime = length(stsTime);
lengthLtsTime = length(ltsTime);

% Frame extraction 
beginningFrame = timingEst;
disp(['Frame beginning (effectively used) :', num2str(beginningFrame)]);
signalCut = signal(beginningFrame:end);
if rc.upsample
    signalCut = downsample(signalCut, rc.USF);
end

lengthFrame = length(ca)...
                +10*lengthStsTime...
                + sc.twoLtsCpLength...
                + 2*lengthLtsTime...
                + sc.signalFieldCpLength...
                + sc.signalFieldLength...
                + sc.numOFDMSymbolsPerFrame*(sc.CPLength+sc.numTotalCarriers);
frame = signalCut(1:lengthFrame);

plotSignalMagnitude(signal, 'Samples', 'Abs val of received signal', beginningFrame, beginningFrame+lengthFrame-1, 'green')
plotSignalMagnitude(frame, 'Samples', 'Received frame')

frameCorrected = frame;

% Carrier frequency offset correction
if rc.cfoCorrection
    % Remove frequency offset and intial phase of the whole frame
    t = (0:length(frameCorrected)-1)*Ts;
    frameCorrected = frameCorrected .* exp(-1j*2*pi*frequencyOffsetEst*t).';
    frameCorrected = frameCorrected * exp(-1j*initialPhaseEst).';
end

if rc.upsample
    figure;
    hax=axes;
    plot(abs(signal))
    title('Abs val of upsampled received signal');
    VL1 = fineTau;
    line([VL1 VL1],get(hax,'YLim'), 'Color', [0 1 0]);
    VL2 = fineTau+rc.USF*lengthFrame-1;
    line([VL2 VL2],get(hax,'YLim'), 'Color', [0 1 0]);
end

plotSignalMagnitude(frameCorrected, 'Samples', 'Corrected received frame')
    
end